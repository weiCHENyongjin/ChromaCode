"""
make_figures.py
===============
Generate the bilingual (English / Chinese) result figures used in the READMEs.

生成 README 使用的中英双语结果图。中文图使用系统中文字体，避免乱码。

Output (in ../figures/):
    fig1_concept_{en,zh}.png        captured signal + its spectrum (5 LED lines)
    fig2_reconstruction_{en,zh}.png 10-channel ground-truth vs reconstructed
    fig3_accuracy_{en,zh}.png       per-channel RMSE and Pearson r

Run:  python code/make_figures.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGDIR = os.path.join(ROOT, "figures")
sys.path.insert(0, HERE)

import iq_sensing_system as sim

# ---- Chinese font (macOS) -------------------------------------------------
_CJK_CANDIDATES = [
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
    "/System/Library/Fonts/STHeiti Medium.ttc",
    "/System/Library/Fonts/PingFang.ttc",
]
ZH_FONT = None
for _p in _CJK_CANDIDATES:
    if os.path.exists(_p):
        fm.fontManager.addfont(_p)
        ZH_FONT = fm.FontProperties(fname=_p).get_name()
        break


def set_lang(lang):
    """Switch matplotlib fonts for the given language."""
    if lang == "zh" and ZH_FONT:
        mpl.rcParams["font.sans-serif"] = [ZH_FONT, "DejaVu Sans"]
    else:
        mpl.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["axes.unicode_minus"] = False


# ---- label dictionaries ---------------------------------------------------
L = {
    "en": {
        "concept_title": "How it works: the spectrum is encoded onto the light source",
        "sig_title": "(a) Captured single-detector signal (first 0.5 s)",
        "sig_x": "time (s)", "sig_y": "intensity (a.u.)",
        "fft_title": "(b) Its spectrum — 5 LED modulation frequencies, each carrying 2 wavelengths",
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


def compute():
    """Run the validated pipeline once; return everything the figures need."""
    t_hi = np.linspace(0, sim.T_SIM, int(sim.T_SIM * sim.OVERFS), endpoint=False)
    true_refl = {wl: sim.target_reflectance(t_hi, wl) for wl in sim.WLS10}
    sig_q, t_sensor, _ = sim.synthesize_sensor_signal(true_refl, seed=42)
    recon, t_valid = sim.reconstruct_all(sig_q, t_sensor)
    metrics = sim.evaluate(recon, t_valid, true_refl, t_hi)
    return dict(t_hi=t_hi, true_refl=true_refl, sig_q=sig_q,
                t_sensor=t_sensor, recon=recon, t_valid=t_valid, metrics=metrics)


# ---- figure builders ------------------------------------------------------
def fig_concept(D, lang):
    set_lang(lang); T = L[lang]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7))

    m = D["t_sensor"] <= 0.5
    ax1.plot(D["t_sensor"][m], D["sig_q"][m], color="#1f77b4", lw=1.2)
    ax1.set_title(T["sig_title"], fontsize=12)
    ax1.set_xlabel(T["sig_x"]); ax1.set_ylabel(T["sig_y"])
    ax1.grid(alpha=0.3)

    sig = D["sig_q"] - D["sig_q"].mean()
    F = np.fft.rfftfreq(sig.size, 1 / sim.FS)
    mag = np.abs(np.fft.rfft(sig)) / sig.size * 2
    band = F <= 70
    ax2.plot(F[band], mag[band], color="#555", lw=1.0)
    for f in sim.FREQS:
        ax2.axvline(f, color="#d62728", ls="--", lw=1, alpha=0.6)
        ax2.annotate(f"{f} Hz", (f, mag[band].max() * 0.92),
                     ha="center", fontsize=9, color="#d62728")
    ax2.set_title(T["fft_title"], fontsize=12)
    ax2.set_xlabel(T["fft_x"]); ax2.set_ylabel(T["fft_y"])
    ax2.grid(alpha=0.3)
    ax2.text(0.99, 0.80, T["fft_note"], transform=ax2.transAxes,
             ha="right", fontsize=9, style="italic", color="#444")

    fig.suptitle(T["concept_title"], fontsize=13, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    _save(fig, f"fig1_concept_{lang}")


def fig_reconstruction(D, lang):
    set_lang(lang); T = L[lang]
    mean_rmse = np.mean([D["metrics"][wl]["rmse"] for wl in sim.WLS10])
    mean_r = np.mean([D["metrics"][wl]["corr"] for wl in sim.WLS10])
    colors = plt.cm.rainbow(np.linspace(0.05, 0.95, 10))

    fig, axes = plt.subplots(5, 2, figsize=(13, 13), sharex=True)
    for k, (f, phi, wl) in enumerate(sim.PAIRS):
        ax = axes[k // 2, k % 2]
        ri = np.interp(D["t_valid"], D["t_hi"], D["true_refl"][wl])
        ax.plot(D["t_valid"], ri, color=colors[k], lw=2.2, alpha=0.6,
                label=T["truth"])
        ax.plot(D["t_valid"], D["recon"][wl], "k--", lw=1.3, label=T["recon"])
        ax.set_ylim(0, 1.05)
        ph = "I" if phi < 0.1 else "Q"
        ax.set_title(f"{wl} nm  |  {f} Hz-{ph}  |  RMSE={D['metrics'][wl]['rmse']:.3f}",
                     fontsize=10)
        ax.grid(alpha=0.25)
        if k == 0:
            ax.legend(fontsize=9, loc="upper right")
        if k % 2 == 0:
            ax.set_ylabel(T["refl_y"], fontsize=9)
    axes[-1, 0].set_xlabel(T["t_x"]); axes[-1, 1].set_xlabel(T["t_x"])
    fig.suptitle(T["recon_title"].format(mean_rmse, mean_r),
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    _save(fig, f"fig2_reconstruction_{lang}")


def fig_accuracy(D, lang):
    set_lang(lang); T = L[lang]
    wls = sim.WLS10
    rmse = [D["metrics"][wl]["rmse"] for wl in wls]
    corr = [D["metrics"][wl]["corr"] for wl in wls]
    x = np.arange(len(wls))
    colors = plt.cm.rainbow(np.linspace(0.05, 0.95, 10))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.6))
    ax1.bar(x, rmse, color=colors)
    ax1.axhline(np.mean(rmse), color="k", ls="--", lw=1,
                label=f"mean = {np.mean(rmse):.3f}")
    ax1.set_xticks(x); ax1.set_xticklabels(wls, rotation=45, fontsize=8)
    ax1.set_title(T["rmse"], fontsize=12); ax1.set_xlabel(T["wl"])
    ax1.legend(fontsize=9); ax1.grid(axis="y", alpha=0.3)

    ax2.bar(x, corr, color=colors)
    ax2.axhline(np.mean(corr), color="k", ls="--", lw=1,
                label=f"mean = {np.mean(corr):.3f}")
    ax2.set_ylim(0.8, 1.0)
    ax2.set_xticks(x); ax2.set_xticklabels(wls, rotation=45, fontsize=8)
    ax2.set_title(T["corr"], fontsize=12); ax2.set_xlabel(T["wl"])
    ax2.legend(fontsize=9); ax2.grid(axis="y", alpha=0.3)

    fig.suptitle(T["acc_title"], fontsize=13, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, f"fig3_accuracy_{lang}")


def _save(fig, name):
    path = os.path.join(FIGDIR, name + ".png")
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved", os.path.relpath(path, ROOT))


if __name__ == "__main__":
    if ZH_FONT:
        print(f"Chinese font: {ZH_FONT}")
    else:
        print("WARNING: no CJK font found; Chinese figures may show boxes.")
    D = compute()
    for lang in ("en", "zh"):
        fig_concept(D, lang)
        fig_reconstruction(D, lang)
        fig_accuracy(D, lang)
    print("done.")
