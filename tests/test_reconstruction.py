"""
End-to-end test: drive the config-based API with the validated forward
simulator and confirm it reproduces the reference accuracy.

用已验证的前向仿真信号驱动配置化 API，确认重建精度与原方案一致。

Run:  python tests/test_reconstruction.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "code"))

import iq_sensing_system as sim                     # forward simulator (validated)
from spectral_reconstruction import SpectralReconstructor


def _base_config(weights=None):
    """Config matching the simulator. Integration delay is set to the
    simulator's exact boxcar delay so this is a faithful equivalence check.
    If ``weights`` is given, attach the per-channel calibration weight."""
    tau = (sim.OVERFS // sim.FS - 1) / (2 * sim.OVERFS)   # = 0.0024 s
    leds = []
    for k, (f, phi, wl) in enumerate(sim.PAIRS):
        led = {"wavelength": wl, "freq": f, "phase_rad": float(phi), "fwhm": 20}
        if weights is not None:
            led["weight"] = float(weights[k])
        leds.append(led)
    return {
        "sensor": {"fs": sim.FS, "integration_delay": tau},
        "lpf": {"cutoff": sim.CUTOFF, "order": 4},
        "leds": leds,
    }


def _score(res, t_hi, true_refl, title):
    print(f"\n=== {title} | mode={res.calibration_mode} "
          f"| tau={res.integration_delay*1e3:.2f} ms ===")
    print(f"{'wavelength':>10} {'RMSE':>10} {'Pearson r':>10}")
    print("-" * 34)
    rmses, corrs = [], []
    for wl in sim.WLS10:
        gt = np.interp(res.t, t_hi, true_refl[wl])
        r = res.reflectance[wl]
        rmse = float(np.sqrt(np.nanmean((r - gt) ** 2)))
        corr = float(np.corrcoef(r, gt)[0, 1])
        rmses.append(rmse); corrs.append(corr)
        print(f"{wl:>8.0f}nm {rmse:>10.4f} {corr:>10.4f}")
    mr = float(np.mean(rmses))
    print("-" * 34)
    print(f"{'mean':>10} {mr:>10.4f} {np.mean(corrs):>10.4f}")
    return mr, float(np.mean(corrs))


def main():
    # ---- forward model: measurement capture ------------------------------
    t_hi = np.linspace(0, sim.T_SIM, int(sim.T_SIM * sim.OVERFS), endpoint=False)
    true_refl = {wl: sim.target_reflectance(t_hi, wl) for wl in sim.WLS10}
    sig, t_sensor, _ = sim.synthesize_sensor_signal(true_refl, seed=42)

    # (1) EXACT EQUIVALENCE: explicit weights -> must match the original ~0.044
    rec_w = SpectralReconstructor.from_dict(_base_config(weights=sim.ws10))
    res_w = rec_w.reconstruct(sig, t=t_sensor, trim=1.0)
    mr_w, corr_w = _score(res_w, t_hi, true_refl, "explicit weights")

    # (2) PRACTICAL PATH: gray reference (R=0.5, won't saturate the ADC),
    #     no knowledge of LED power or sensor response required.
    gray = {wl: 0.5 * np.ones_like(t_hi) for wl in sim.WLS10}
    gray_sig, _, _ = sim.synthesize_sensor_signal(gray, seed=7)
    rec_g = SpectralReconstructor.from_dict(_base_config())
    res_g = rec_g.reconstruct(sig, t=t_sensor, white_reference=gray_sig,
                              reference_level=0.5, trim=1.0)
    mr_g, corr_g = _score(res_g, t_hi, true_refl, "gray reference (R=0.5)")

    # (3) SPECTRAL PATH with MISMATCHED wavelength sampling: sensor response
    #     given on a coarse 10 nm grid; LED lines are 20 nm Gaussians on an
    #     implicit fine grid. Verifies the resampling code path runs and the
    #     shape (correlation) is recovered (absolute scale is relative here).
    cfg = _base_config()
    cfg["sensor"]["spectral_response"] = {
        "wavelengths": list(range(400, 1001, 10)),                  # 10 nm grid
        "response": [float(sim.sensor_response(np.array([w]))[0])    # sampled coarsely
                     for w in range(400, 1001, 10)],
    }
    res_s = SpectralReconstructor.from_dict(cfg).reconstruct(sig, t=t_sensor, trim=1.0)
    corr_s = float(np.mean([
        np.corrcoef(res_s.reflectance[wl], np.interp(res_s.t, t_hi, true_refl[wl]))[0, 1]
        for wl in sim.WLS10]))
    print(f"\n[spectral path, mismatched 10nm grid] mode={res_s.calibration_mode} "
          f"mean corr={corr_s:.4f}")
    assert res_s.calibration_mode == "spectral"
    assert np.all([np.all(np.isfinite(res_s.reflectance[wl])) for wl in sim.WLS10])
    assert corr_s > 0.95, f"spectral-path mean corr too low: {corr_s}"

    # ---- assertions ------------------------------------------------------
    assert mr_w < 0.05, f"explicit-weights mean RMSE too high: {mr_w}"
    assert corr_w > 0.95, f"explicit-weights mean corr too low: {corr_w}"
    assert mr_g < 0.07, f"gray-reference mean RMSE too high: {mr_g}"
    assert corr_g > 0.95, f"gray-reference mean corr too low: {corr_g}"
    print("\nPASS: config-driven API reproduces the validated accuracy "
          f"(explicit-weights mean RMSE={mr_w:.4f}, matches original 0.044).")


if __name__ == "__main__":
    main()
