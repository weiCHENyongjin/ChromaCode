"""Generate the bilingual (English / Chinese) result figures for the READMEs.

生成 README 使用的中英双语结果图。中文图使用系统中文字体，避免乱码。

Output (in ../figures/)
-----------------------
``fig0_spectra_{en,zh}.png``
    LED emission bands + sensor spectral response.
``fig1_concept_{en,zh}.png``
    Captured signal + its spectrum (the five LED modulation lines).
``fig2_reconstruction_{en,zh}.png``
    10-channel ground-truth vs reconstructed reflectance.
``fig3_accuracy_{en,zh}.png``
    Per-channel RMSE and Pearson r.

Run ``python code/make_figures.py``. See ``../CODE_STYLE.md`` for conventions.
"""

from __future__ import annotations

import os
import sys
from typing import Dict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGURE_DIR = os.path.join(ROOT, "figures")
sys.path.insert(0, HERE)

from iq_sensing_system import (  # noqa: E402  (local import after path setup)
    MODULATION_FREQUENCIES_HZ,
    SPECTRAL_GRID_NM,
    WAVELENGTHS_NM,
    ForwardSimulator,
    evaluate_reconstruction,
    synthetic_target_reflectance,
)
from spectral_reconstruction import (  # noqa: E402
    gaussian_emission,
    silicon_sensor_response,
)


# ─────────────────────────────────────────────────────────────────────────
# Fonts
# ─────────────────────────────────────────────────────────────────────────
_CJK_FONT_CANDIDATES = [
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
    "/System/Library/Fonts/STHeiti Medium.ttc",
    "/System/Library/Fonts/PingFang.ttc",
]
CJK_FONT_NAME = None
for _path in _CJK_FONT_CANDIDATES:
    if os.path.exists(_path):
        fm.fontManager.addfont(_path)
        CJK_FONT_NAME = fm.FontProperties(fname=_path).get_name()
        break


def set_language(lang: str) -> None:
    """Switch matplotlib fonts for the given language ('en' or 'zh')."""
    if lang == "zh" and CJK_FONT_NAME:
        mpl.rcParams["font.sans-serif"] = [CJK_FONT_NAME, "DejaVu Sans"]
    else:
        mpl.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["axes.unicode_minus"] = False


# ─────────────────────────────────────────────────────────────────────────
# Label dictionaries
# ─────────────────────────────────────────────────────────────────────────
LABELS = {
    "en": {
        "spectra_title": "LED emission bands and sensor spectral response "
                         "(10 channels = 5 frequencies x IQ phases, d=0.5)",
        "spectra_x": "wavelength (nm)", "spectra_y": "normalized intensity / response",
        "sensor_resp": "sensor response",
        "concept_title": "How it works: the spectrum is encoded onto the light source",
        "sig_title": "(a) Captured single-detector signal (first 0.5 s)",
        "sig_x": "time (s)", "sig_y": "intensity (a.u.)",
        "fft_title": "(b) Its spectrum - 5 LED modulation frequencies, each carrying 2 wavelengths",
        "fft_x": "frequency (Hz)", "fft_y": "magnitude",
        "fft_note": "lock-in demodulation separates each line",
        "recon_title": "10-channel reflectance reconstruction (mean RMSE = {:.3f}, r = {:.3f})",
        "truth": "ground truth", "recon": "reconstructed",
        "t_x": "time (s)", "refl_y": "reflectance",
        "acc_title": "Per-channel accuracy",
        "rmse": "RMSE (lower is better)", "corr": "Pearson r (higher is better)",
        "wl": "wavelength (nm)",
    },
    "zh": {
        "spectra_title": "LED 发射谱与传感器光谱响应（10 通道 = 5 频率 x IQ 相位，d=0.5）",
        "spectra_x": "波长 (nm)", "spectra_y": "归一化强度 / 响应",
        "sensor_resp": "传感器响应",
        "concept_title": "原理示意：光谱信息被编码到光源上",
        "sig_title": "(a) 单传感器采集信号（前 0.5 秒）",
        "sig_x": "时间 (s)", "sig_y": "光强 (a.u.)",
        "fft_title": "(b) 其频谱 —— 5 个 LED 调制频率，每个承载 2 个波长",
        "fft_x": "频率 (Hz)", "fft_y": "幅值",
        "fft_note": "锁相解调逐一分离各谱线",
        "recon_title": "10 通道反射率重建（均值 RMSE = {:.3f}，r = {:.3f}）",
        "truth": "真实值", "recon": "重建值",
        "t_x": "时间 (s)", "refl_y": "反射率",
        "acc_title": "各通道精度",
        "rmse": "RMSE（越低越好）", "corr": "皮尔逊 r（越高越好）",
        "wl": "波长 (nm)",
    },
}


# ─────────────────────────────────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────────────────────────────────
def compute() -> dict:
    """Run the validated pipeline once; return everything the figures need."""
    simulator = ForwardSimulator.default()
    time_hi = simulator.oversample_time_s()
    true_reflectance = {wl: synthetic_target_reflectance(time_hi, wl)
                        for wl in WAVELENGTHS_NM}
    signal_q = simulator.synthesize(true_reflectance)
    sensor_time = simulator.sensor_time_s()
    result = simulator.reconstructor().reconstruct(
        signal_q, time_s=sensor_time, trim_s=1.0)
    metrics = evaluate_reconstruction(result, time_hi, true_reflectance)
    return dict(simulator=simulator, time_hi=time_hi,
                true_reflectance=true_reflectance, signal_q=signal_q,
                sensor_time=sensor_time, result=result, metrics=metrics)


# ─────────────────────────────────────────────────────────────────────────
# Figure builders
# ─────────────────────────────────────────────────────────────────────────
def _channel_color(index: int) -> np.ndarray:
    return plt.cm.rainbow(np.linspace(0.05, 0.95, 10))[index]


def figure_spectra(data: dict, lang: str) -> None:
    set_language(lang)
    text = LABELS[lang]
    simulator = data["simulator"]
    response = silicon_sensor_response(SPECTRAL_GRID_NM)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(SPECTRAL_GRID_NM, response, "k-", lw=2, label=text["sensor_resp"])
    for index, channel in enumerate(simulator.channels):
        band = gaussian_emission(SPECTRAL_GRID_NM, channel.wavelength_nm,
                                 channel.fwhm_nm)
        color = _channel_color(index)
        ax.fill_between(SPECTRAL_GRID_NM, band, color=color, alpha=0.45)
        ax.plot(SPECTRAL_GRID_NM, band, color=color, lw=1)
        phase = "I" if channel.phase_rad < 0.1 else "Q"
        ax.annotate(f"{channel.wavelength_nm:.0f}nm\n{channel.frequency_hz:.0f}Hz-{phase}",
                    (channel.wavelength_nm, 1.04), ha="center", va="bottom",
                    fontsize=7.5, color=color)
    ax.set_xlim(400, 1000)
    ax.set_ylim(0, 1.22)
    ax.set_xlabel(text["spectra_x"])
    ax.set_ylabel(text["spectra_y"])
    ax.set_title(text["spectra_title"], fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    _save(fig, f"fig0_spectra_{lang}")


def figure_concept(data: dict, lang: str) -> None:
    set_language(lang)
    text = LABELS[lang]
    sample_rate_hz = data["simulator"].params.sample_rate_hz
    sensor_time = data["sensor_time"]
    signal_q = data["signal_q"]

    fig, (ax_signal, ax_spectrum) = plt.subplots(2, 1, figsize=(11, 7))

    window = sensor_time <= 0.5
    ax_signal.plot(sensor_time[window], signal_q[window], color="#1f77b4", lw=1.2)
    ax_signal.set_title(text["sig_title"], fontsize=12)
    ax_signal.set_xlabel(text["sig_x"])
    ax_signal.set_ylabel(text["sig_y"])
    ax_signal.grid(alpha=0.3)

    centered = signal_q - signal_q.mean()
    freqs = np.fft.rfftfreq(centered.size, 1 / sample_rate_hz)
    magnitude = np.abs(np.fft.rfft(centered)) / centered.size * 2
    band = freqs <= 70
    ax_spectrum.plot(freqs[band], magnitude[band], color="#555", lw=1.0)
    for frequency_hz in MODULATION_FREQUENCIES_HZ:
        ax_spectrum.axvline(frequency_hz, color="#d62728", ls="--", lw=1, alpha=0.6)
        ax_spectrum.annotate(f"{frequency_hz:.0f} Hz",
                             (frequency_hz, magnitude[band].max() * 0.92),
                             ha="center", fontsize=9, color="#d62728")
    ax_spectrum.set_title(text["fft_title"], fontsize=12)
    ax_spectrum.set_xlabel(text["fft_x"])
    ax_spectrum.set_ylabel(text["fft_y"])
    ax_spectrum.grid(alpha=0.3)
    ax_spectrum.text(0.99, 0.80, text["fft_note"], transform=ax_spectrum.transAxes,
                     ha="right", fontsize=9, style="italic", color="#444")

    fig.suptitle(text["concept_title"], fontsize=13, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    _save(fig, f"fig1_concept_{lang}")


def figure_reconstruction(data: dict, lang: str) -> None:
    set_language(lang)
    text = LABELS[lang]
    simulator, result, metrics = data["simulator"], data["result"], data["metrics"]
    time_hi, true_reflectance = data["time_hi"], data["true_reflectance"]

    mean_rmse = np.mean([m["rmse"] for m in metrics.values()])
    mean_corr = np.mean([m["corr"] for m in metrics.values()])

    fig, axes = plt.subplots(5, 2, figsize=(13, 13), sharex=True)
    for index, channel in enumerate(simulator.channels):
        ax = axes[index // 2, index % 2]
        wavelength_nm = channel.wavelength_nm
        truth = np.interp(result.time_s, time_hi, true_reflectance[wavelength_nm])
        ax.plot(result.time_s, truth, color=_channel_color(index), lw=2.2,
                alpha=0.6, label=text["truth"])
        ax.plot(result.time_s, result.reflectance[wavelength_nm], "k--", lw=1.3,
                label=text["recon"])
        ax.set_ylim(0, 1.05)
        phase = "I" if channel.phase_rad < 0.1 else "Q"
        ax.set_title(f"{wavelength_nm:.0f} nm  |  {channel.frequency_hz:.0f} Hz-{phase}"
                     f"  |  RMSE={metrics[wavelength_nm]['rmse']:.3f}", fontsize=10)
        ax.grid(alpha=0.25)
        if index == 0:
            ax.legend(fontsize=9, loc="upper right")
        if index % 2 == 0:
            ax.set_ylabel(text["refl_y"], fontsize=9)
    axes[-1, 0].set_xlabel(text["t_x"])
    axes[-1, 1].set_xlabel(text["t_x"])
    fig.suptitle(text["recon_title"].format(mean_rmse, mean_corr),
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    _save(fig, f"fig2_reconstruction_{lang}")


def figure_accuracy(data: dict, lang: str) -> None:
    set_language(lang)
    text = LABELS[lang]
    metrics = data["metrics"]
    wavelengths = WAVELENGTHS_NM
    rmse = [metrics[wl]["rmse"] for wl in wavelengths]
    corr = [metrics[wl]["corr"] for wl in wavelengths]
    positions = np.arange(len(wavelengths))
    colors = plt.cm.rainbow(np.linspace(0.05, 0.95, 10))

    fig, (ax_rmse, ax_corr) = plt.subplots(1, 2, figsize=(13, 4.6))
    ax_rmse.bar(positions, rmse, color=colors)
    ax_rmse.axhline(np.mean(rmse), color="k", ls="--", lw=1,
                    label=f"mean = {np.mean(rmse):.3f}")
    ax_rmse.set_xticks(positions)
    ax_rmse.set_xticklabels([f"{wl:.0f}" for wl in wavelengths], rotation=45, fontsize=8)
    ax_rmse.set_title(text["rmse"], fontsize=12)
    ax_rmse.set_xlabel(text["wl"])
    ax_rmse.legend(fontsize=9)
    ax_rmse.grid(axis="y", alpha=0.3)

    ax_corr.bar(positions, corr, color=colors)
    ax_corr.axhline(np.mean(corr), color="k", ls="--", lw=1,
                    label=f"mean = {np.mean(corr):.3f}")
    ax_corr.set_ylim(0.8, 1.0)
    ax_corr.set_xticks(positions)
    ax_corr.set_xticklabels([f"{wl:.0f}" for wl in wavelengths], rotation=45, fontsize=8)
    ax_corr.set_title(text["corr"], fontsize=12)
    ax_corr.set_xlabel(text["wl"])
    ax_corr.legend(fontsize=9)
    ax_corr.grid(axis="y", alpha=0.3)

    fig.suptitle(text["acc_title"], fontsize=13, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, f"fig3_accuracy_{lang}")


def _save(fig, name: str) -> None:
    path = os.path.join(FIGURE_DIR, name + ".png")
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved", os.path.relpath(path, ROOT))


def main() -> None:
    if CJK_FONT_NAME:
        print(f"Chinese font: {CJK_FONT_NAME}")
    else:
        print("WARNING: no CJK font found; Chinese figures may show boxes.")
    data = compute()
    for lang in ("en", "zh"):
        figure_spectra(data, lang)
        figure_concept(data, lang)
        figure_reconstruction(data, lang)
        figure_accuracy(data, lang)
    print("done.")


if __name__ == "__main__":
    main()
