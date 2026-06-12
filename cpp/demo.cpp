// ChromaCode C++ demo / validation.
//
// Loads the SAME JSON config the Python reconstructor uses, plus a captured
// signal, reconstructs the wavelength channels, and scores each channel's RMSE
// and Pearson r against the analytic ground truth.
//
// Build (see cpp/README.md):
//   clang++ -std=c++17 -O2 -I$(brew --prefix armadillo)/include \
//     demo.cpp -L$(brew --prefix armadillo)/lib -larmadillo -o chromacode_demo
//
// Run:
//   ./chromacode_demo ../config/default_10ch.json sample_signal.csv

#include <cmath>
#include <cstdio>

#include "chromacode.hpp"

using chromacode::ReconstructionConfig;
using chromacode::SpectralReconstructor;
using chromacode::load_config_file;

namespace {

// Analytic ground-truth reflectance (matches synthetic_target_reflectance).
arma::vec ground_truth(const arma::vec& t, double wl) {
    const double f1 = 0.5 + (wl - 450.0) / 400.0;
    const double f2 = 1.7 * f1;
    arma::vec r = 0.4 + 0.3 * std::exp(-std::pow((wl - 600.0) / 150.0, 2))
                + 0.15 * arma::sin(2.0 * M_PI * f1 * t)
                + 0.08 * arma::sin(2.0 * M_PI * f2 * t + M_PI / 3.0);
    return arma::clamp(r, 0.05, 0.95);
}

double pearson(const arma::vec& a, const arma::vec& b) {
    return arma::as_scalar(arma::cor(a, b));
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
                     "usage: %s <config.json> <signal.csv> [white_reference.csv]\n",
                     argv[0]);
        return 1;
    }

    ReconstructionConfig config;
    try {
        config = load_config_file(argv[1]);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error loading config: %s\n", e.what());
        return 1;
    }

    arma::vec signal;
    if (!signal.load(argv[2], arma::raw_ascii) || signal.is_empty()) {
        std::fprintf(stderr, "error: cannot load signal from %s\n", argv[2]);
        return 1;
    }
    arma::vec white;
    if (argc >= 4 && !white.load(argv[3], arma::raw_ascii)) {
        std::fprintf(stderr, "error: cannot load white reference from %s\n", argv[3]);
        return 1;
    }

    SpectralReconstructor reconstructor(config);
    auto result = reconstructor.reconstruct(signal, white,
                                            /*reference_level=*/1.0,
                                            /*trim_s=*/1.0);

    std::printf("config           : %s\n", argv[1]);
    std::printf("calibration mode : %s\n", result.calibration_mode.c_str());
    std::printf("integration delay: %.3f ms\n", result.integration_delay_s * 1e3);
    std::printf("%10s %10s %10s\n", "wavelength", "RMSE", "Pearson r");
    std::printf("----------------------------------\n");

    double sum_rmse = 0.0, sum_corr = 0.0;
    int n = 0;
    for (double wl : result.wavelengths_nm) {
        const arma::vec& R = result.reflectance[wl];
        arma::vec gt = ground_truth(result.time_s, wl);
        double rmse = std::sqrt(arma::mean(arma::square(R - gt)));
        double corr = pearson(R, gt);
        std::printf("%8.0fnm %10.4f %10.4f\n", wl, rmse, corr);
        sum_rmse += rmse;
        sum_corr += corr;
        ++n;
    }
    const double mean_rmse = sum_rmse / n;
    const double mean_corr = sum_corr / n;
    std::printf("----------------------------------\n");
    std::printf("%10s %10.4f %10.4f\n", "mean", mean_rmse, mean_corr);

    if (mean_rmse < 0.05 && mean_corr > 0.95) {
        std::printf("\nPASS: C++ port reproduces the validated accuracy "
                    "(mean RMSE %.4f ~ 0.044).\n", mean_rmse);
        return 0;
    }
    std::printf("\nFAIL: mean RMSE %.4f / corr %.4f out of expected range.\n",
                mean_rmse, mean_corr);
    return 2;
}
