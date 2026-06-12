"""Minimal end-to-end example of the config-driven reconstruction API.

最小可运行示例：从 YAML 配置文件加载参数，传入采集信号，重建各波长反射率。

The validated simulator is used only to *produce a demo capture*; in real use
replace :func:`demo_capture` with your detector signal (e.g. loaded from a
``.npy`` / ``.csv`` file).

Run:  python examples/example_usage.py
"""

from __future__ import annotations

import os
import sys
from typing import Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))   # python/examples
PYTHON_ROOT = os.path.dirname(HERE)                  # python/
REPO_ROOT = os.path.dirname(PYTHON_ROOT)             # repo root
sys.path.insert(0, PYTHON_ROOT)

from iq_sensing_system import (  # noqa: E402
    WAVELENGTHS_NM,
    ForwardSimulator,
    synthetic_target_reflectance,
)
from spectral_reconstruction import SpectralReconstructor  # noqa: E402


def demo_capture() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    """Return (signal, gray_reference, sensor_time, time_hi, ground_truth)."""
    simulator = ForwardSimulator.default()
    time_hi = simulator.oversample_time_s()
    sensor_time = simulator.sensor_time_s()
    true_reflectance = {wl: synthetic_target_reflectance(time_hi, wl)
                        for wl in WAVELENGTHS_NM}
    gray = {wl: 0.5 * np.ones_like(time_hi) for wl in WAVELENGTHS_NM}  # R=0.5 card
    signal = simulator.synthesize(true_reflectance, seed=42)
    gray_signal = simulator.synthesize(gray, seed=7)
    return signal, gray_signal, sensor_time, time_hi, true_reflectance


def main() -> None:
    signal, gray_signal, sensor_time, time_hi, true_reflectance = demo_capture()

    # 1) Load the system description from a config file.
    config_path = os.path.join(REPO_ROOT, "config", "example_config.yaml")
    reconstructor = SpectralReconstructor.from_config_file(config_path)

    # 2) Reconstruct. A gray reference (known reflectance 0.5) gives absolute
    #    calibration without needing the sensor spectral response.
    result = reconstructor.reconstruct(
        signal, time_s=sensor_time, white_reference=gray_signal,
        reference_level=0.5, trim_s=1.0)

    # 3) Use the results.
    print(f"calibration mode : {result.calibration_mode}")
    print(f"integration delay: {result.integration_delay_s * 1e3:.2f} ms")
    print(f"reconstructed {result.meta['n_channels']} wavelength channels:\n")
    print(f"{'wavelength':>10} {'mean R':>8} {'RMSE':>8}")
    for wavelength_nm in result.wavelengths_nm:
        series = result.reflectance[wavelength_nm]
        truth = np.interp(result.time_s, time_hi, true_reflectance[wavelength_nm])
        rmse = np.sqrt(np.nanmean((series - truth) ** 2))
        print(f"{wavelength_nm:>8.0f}nm {np.nanmean(series):>8.3f} {rmse:>8.4f}")

    # Optional: plot if matplotlib is available.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        wavelengths, _ = result.as_matrix()
        fig, axes = plt.subplots(len(wavelengths), 1, figsize=(11, 16), sharex=True)
        for ax, wavelength_nm in zip(axes, wavelengths):
            truth = np.interp(result.time_s, time_hi, true_reflectance[wavelength_nm])
            ax.plot(result.time_s, truth, lw=2, alpha=0.5, label="ground truth")
            ax.plot(result.time_s, result.reflectance[wavelength_nm], "k--", lw=1.2,
                    label="reconstructed")
            ax.set_ylabel(f"{wavelength_nm:.0f} nm", fontsize=8)
            ax.set_ylim(0, 1.1)
        axes[0].legend(fontsize=8, loc="upper right")
        axes[-1].set_xlabel("time (s)")
        out_path = os.path.join(HERE, "example_reconstruction.png")
        fig.suptitle("Config-driven reconstruction (gray-reference calibration)")
        plt.tight_layout()
        plt.savefig(out_path, dpi=130)
        plt.close(fig)
        print(f"\nfigure saved: {out_path}")
    except Exception as exc:  # pragma: no cover - plotting is optional
        print(f"(plot skipped: {exc})")


if __name__ == "__main__":
    main()
