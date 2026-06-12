"""End-to-end test of the config-driven reconstruction API.

用已验证的前向仿真信号驱动配置化 API，确认重建精度与原方案一致。

Three calibration paths are exercised:
    1. explicit weights  -> exact reproduction of the published ~0.044 RMSE
    2. gray reference    -> practical, no LED-power / response knowledge
    3. spectral response on a *mismatched* wavelength grid -> resampling path

Run:  python tests/test_reconstruction.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))   # python/tests
PYTHON_ROOT = os.path.dirname(HERE)                  # python/
sys.path.insert(0, PYTHON_ROOT)

from iq_sensing_system import (  # noqa: E402
    SPECTRAL_GRID_NM,
    WAVELENGTHS_NM,
    ForwardSimulator,
    evaluate_reconstruction,
    synthetic_target_reflectance,
)
from spectral_reconstruction import (  # noqa: E402
    ReconstructionConfig,
    SpectralReconstructor,
    silicon_sensor_response,
)


def _score(result, time_hi, true_reflectance, title):
    """Print and return (mean RMSE, mean Pearson r) for a reconstruction."""
    metrics = evaluate_reconstruction(result, time_hi, true_reflectance)
    print(f"\n=== {title} | mode={result.calibration_mode} "
          f"| tau={result.integration_delay_s * 1e3:.2f} ms ===")
    print(f"{'wavelength':>10} {'RMSE':>10} {'Pearson r':>10}")
    print("-" * 34)
    for wl in WAVELENGTHS_NM:
        print(f"{wl:>8.0f}nm {metrics[wl]['rmse']:>10.4f} {metrics[wl]['corr']:>10.4f}")
    mean_rmse = float(np.mean([metrics[wl]["rmse"] for wl in WAVELENGTHS_NM]))
    mean_corr = float(np.mean([metrics[wl]["corr"] for wl in WAVELENGTHS_NM]))
    print("-" * 34)
    print(f"{'mean':>10} {mean_rmse:>10.4f} {mean_corr:>10.4f}")
    return mean_rmse, mean_corr


def main() -> None:
    simulator = ForwardSimulator.default()
    time_hi = simulator.oversample_time_s()
    sensor_time = simulator.sensor_time_s()
    true_reflectance = {wl: synthetic_target_reflectance(time_hi, wl)
                        for wl in WAVELENGTHS_NM}
    signal_q = simulator.synthesize(true_reflectance, seed=42)

    # (1) EXACT EQUIVALENCE: the simulator's matched reconstructor uses the
    #     known per-channel weights -> reproduces the published ~0.044 RMSE.
    result_w = simulator.reconstructor().reconstruct(
        signal_q, time_s=sensor_time, trim_s=1.0)
    mean_rmse_w, mean_corr_w = _score(result_w, time_hi, true_reflectance,
                                      "explicit weights")

    # (2) PRACTICAL PATH: a gray reference (R=0.5, won't saturate the ADC); no
    #     knowledge of LED power or sensor response is required.
    gray = {wl: 0.5 * np.ones_like(time_hi) for wl in WAVELENGTHS_NM}
    gray_signal = simulator.synthesize(gray, seed=7)
    reconstructor = SpectralReconstructor(simulator.reconstruction_config())
    result_g = reconstructor.reconstruct(
        signal_q, time_s=sensor_time, white_reference=gray_signal,
        reference_level=0.5, trim_s=1.0)
    mean_rmse_g, mean_corr_g = _score(result_g, time_hi, true_reflectance,
                                      "gray reference (R=0.5)")

    # (3) SPECTRAL PATH with MISMATCHED wavelength sampling: the sensor response
    #     is supplied on a coarse 10 nm grid while LED lines live on a fine
    #     implicit grid. Verifies the resampling path runs and recovers shape.
    coarse_grid = np.arange(400, 1001, 10)
    config_dict = {
        "sensor": {
            "fs": simulator.params.sample_rate_hz,
            "integration_delay": simulator.params.boxcar_delay_s,
            "spectral_response": {
                "wavelengths": [float(w) for w in coarse_grid],
                "response": [float(silicon_sensor_response(np.array([w]))[0])
                             for w in coarse_grid],
            },
        },
        "lpf": {"cutoff": simulator.params.lpf_cutoff_hz, "order": 4},
        "leds": [{"wavelength": ch.wavelength_nm, "freq": ch.frequency_hz,
                  "phase_rad": ch.phase_rad, "fwhm": ch.fwhm_nm}
                 for ch in simulator.channels],
    }
    result_s = SpectralReconstructor.from_dict(config_dict).reconstruct(
        signal_q, time_s=sensor_time, trim_s=1.0)
    metrics_s = evaluate_reconstruction(result_s, time_hi, true_reflectance)
    mean_corr_s = float(np.mean([metrics_s[wl]["corr"] for wl in WAVELENGTHS_NM]))
    print(f"\n[spectral path, mismatched 10nm grid] mode={result_s.calibration_mode} "
          f"mean corr={mean_corr_s:.4f}")

    # ---- assertions ------------------------------------------------------
    assert isinstance(simulator.reconstruction_config(), ReconstructionConfig)
    assert mean_rmse_w < 0.05, f"explicit-weights mean RMSE too high: {mean_rmse_w}"
    assert mean_corr_w > 0.95, f"explicit-weights mean corr too low: {mean_corr_w}"
    assert mean_rmse_g < 0.07, f"gray-reference mean RMSE too high: {mean_rmse_g}"
    assert mean_corr_g > 0.95, f"gray-reference mean corr too low: {mean_corr_g}"
    assert result_s.calibration_mode == "spectral"
    assert all(np.all(np.isfinite(result_s.reflectance[wl])) for wl in WAVELENGTHS_NM)
    assert mean_corr_s > 0.95, f"spectral-path mean corr too low: {mean_corr_s}"

    print("\nPASS: config-driven API reproduces the validated accuracy "
          f"(explicit-weights mean RMSE={mean_rmse_w:.4f}, matches original 0.044).")


if __name__ == "__main__":
    main()
