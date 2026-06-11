"""
example_usage.py
================
Minimal end-to-end example of the config-driven reconstruction API.

最小可运行示例：从 YAML 配置文件加载参数，传入采集信号，重建各波长反射率。

It uses the validated simulator only to *produce a demo signal*; in real use
you would replace ``demo_signal()`` with your detector capture (e.g. loaded
from a .npy/.csv file).

Run:  python examples/example_usage.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "code"))

from spectral_reconstruction import SpectralReconstructor
import iq_sensing_system as sim   # only to fabricate a demo capture


def demo_signal():
    """Return (signal, gray_reference, t, ground_truth) for the demo."""
    t_hi = np.linspace(0, sim.T_SIM, int(sim.T_SIM * sim.OVERFS), endpoint=False)
    true_refl = {wl: sim.target_reflectance(t_hi, wl) for wl in sim.WLS10}
    gray = {wl: 0.5 * np.ones_like(t_hi) for wl in sim.WLS10}   # R=0.5 gray card
    sig, t_sensor, _ = sim.synthesize_sensor_signal(true_refl, seed=42)
    gray_sig, _, _ = sim.synthesize_sensor_signal(gray, seed=7)
    return sig, gray_sig, t_sensor, t_hi, true_refl


def main():
    sig, gray_sig, t_sensor, t_hi, true_refl = demo_signal()

    # 1) Load the system description from a config file.
    cfg_path = os.path.join(ROOT, "config", "example_config.yaml")
    rec = SpectralReconstructor.from_config_file(cfg_path)

    # 2) Reconstruct. Pass a gray reference (known reflectance 0.5) for
    #    absolute calibration without needing the sensor spectral response.
    res = rec.reconstruct(sig, t=t_sensor,
                          white_reference=gray_sig, reference_level=0.5,
                          trim=1.0)

    # 3) Use the results.
    print(f"calibration mode : {res.calibration_mode}")
    print(f"integration delay: {res.integration_delay*1e3:.2f} ms")
    print(f"reconstructed {res.meta['n_channels']} wavelength channels:\n")
    print(f"{'wavelength':>10} {'mean R':>8} {'RMSE':>8}")
    for wl in res.wavelengths:
        r = res.reflectance[wl]
        gt = np.interp(res.t, t_hi, true_refl[wl])
        rmse = np.sqrt(np.nanmean((r - gt) ** 2))
        print(f"{wl:>8.0f}nm {np.nanmean(r):>8.3f} {rmse:>8.4f}")

    # Optional: plot if matplotlib is available.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        wls, M = res.as_matrix()
        fig, axes = plt.subplots(len(wls), 1, figsize=(11, 16), sharex=True)
        for ax, wl in zip(axes, wls):
            gt = np.interp(res.t, t_hi, true_refl[wl])
            ax.plot(res.t, gt, lw=2, alpha=0.5, label="ground truth")
            ax.plot(res.t, res.reflectance[wl], "k--", lw=1.2, label="reconstructed")
            ax.set_ylabel(f"{wl:.0f} nm", fontsize=8)
            ax.set_ylim(0, 1.1)
        axes[0].legend(fontsize=8, loc="upper right")
        axes[-1].set_xlabel("time (s)")
        out = os.path.join(HERE, "example_reconstruction.png")
        fig.suptitle("Config-driven reconstruction (gray-reference calibration)")
        plt.tight_layout()
        plt.savefig(out, dpi=130)
        print(f"\nfigure saved: {out}")
    except Exception as e:  # pragma: no cover
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
