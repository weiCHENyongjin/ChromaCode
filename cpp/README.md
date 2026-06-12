# ChromaCode — C++ reconstruction core

A self-contained C++17 port of the Python `SpectralReconstructor` inverse path:
**reference → mix → zero-phase Butterworth low-pass → calibrate**. It reproduces
the Python reference accuracy channel-for-channel (mean RMSE ≈ 0.044).

## Files

| File | Purpose |
|------|---------|
| `chromacode.hpp` | Header-only library: Butterworth low-pass design, zero-phase `sosfiltfilt`, `SpectralReconstructor`. |
| `demo.cpp` | Loads a captured signal, reconstructs the 10 channels, scores RMSE / Pearson r against the analytic ground truth. |
| `sample_signal.csv` | A captured signal (2400 samples) exported from the Python simulator (`seed=42`). |
| `CMakeLists.txt` | Optional CMake build. |

## Dependencies

- **[Armadillo](https://arma.sourceforge.net/)** (+ BLAS/LAPACK) — vectors and
  element-wise math. The Butterworth design and zero-phase filter are
  implemented in `chromacode.hpp` (no DSP library needed).

```bash
brew install armadillo          # macOS
# sudo apt install libarmadillo-dev   # Debian/Ubuntu
```

> `dlib` / `NLopt` / `Eigen` are not required by the core; they are natural fits
> for future extensions (e.g. optimization-based calibration or multi-pixel work).

## Build & run

**Direct (clang++ / g++):**

```bash
ARMA=$(brew --prefix armadillo)
clang++ -std=c++17 -O2 -I"$ARMA/include" demo.cpp \
        -L"$ARMA/lib" -larmadillo -o chromacode_demo
./chromacode_demo sample_signal.csv
```

**CMake:**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/chromacode_demo sample_signal.csv
```

Expected output (matches the Python reference exactly):

```
calibration mode : weights
integration delay: 2.400 ms
wavelength       RMSE  Pearson r
----------------------------------
     450nm     0.0246     0.9791
     ...
     850nm     0.0717     0.9807
----------------------------------
      mean     0.0439     0.9642

PASS: C++ port reproduces the validated accuracy (mean RMSE 0.0439 ~ 0.044).
```

## Library usage

```cpp
#include "chromacode.hpp"
using namespace chromacode;

ReconstructionConfig cfg;
cfg.sample_rate_hz      = 200.0;
cfg.lpf_cutoff_hz       = 3.5;
cfg.lpf_order           = 4;
cfg.integration_delay_s = 0.0024;          // or set integration_time_s
cfg.channels = {                            // {wavelength, freq, phase, weight}
    {450, 13, 0.0,        0.1167},
    {494, 13, M_PI / 2.0, 0.1880},
    // ... 10 channels
};

SpectralReconstructor rec(cfg);
arma::vec signal;  signal.load("signal.csv", arma::raw_ascii);

// Weights mode (per-channel weight known):
auto result = rec.reconstruct(signal, /*white=*/arma::vec(), 1.0, /*trim_s=*/1.0);

// White/gray-reference mode (absolute reflectance, response unknown):
//   auto result = rec.reconstruct(signal, white_capture, /*reference_level=*/0.5, 1.0);

const arma::vec& r850 = result.reflectance[850.0];
```

## Validation

`demo.cpp` doubles as a regression test: it asserts `mean RMSE < 0.05` and
`mean r > 0.95`, and prints `PASS`/`FAIL`. The numbers are identical to the
Python `tests/test_reconstruction.py` reference, confirming the port is exact.

The shared theory (definitions, lemmas, proofs) is in
[`../docs/mathematical_theory.md`](../docs/mathematical_theory.md).
