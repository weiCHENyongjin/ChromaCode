"""
多波长 LED IQ 相位编码光强采集与反射率重建系统
Multi-Wavelength LED IQ Phase-Encoded Reflectance Sensing System

最优方案配置（v1）：
  - 传感器：Fs=200Hz，8-bit，SNR≈42dB
  - LED：[13,23,33,43,53]Hz × d=0.5 × 相位{0°,90°} = 10通道
  - 波长：[450,494,539,583,628,672,717,761,806,850] nm
  - 解调：4阶Butterworth LPF（截止3.5Hz）+ 相位时延补偿
  - 性能：均值RMSE=0.044，皮尔逊r=0.964

依赖：numpy, scipy, matplotlib
"""

import numpy as np
import scipy.signal as signal
from scipy.signal import butter, sosfiltfilt
import matplotlib.pyplot as plt
import os

# ─────────────────────────────────────────────────────
# 1. 系统参数（最优配置）
# ─────────────────────────────────────────────────────

FS        = 200          # 传感器采样率 (Hz)
OVERFS    = 5000         # 过采样仿真率 (Hz)  —— 仿真专用，真实系统无需
T_SIM     = 12.0         # 仿真时长 (s)
CUTOFF    = 3.5          # 锁相低通截止频率 (Hz)
SNR_DB    = 42.0         # 信噪比 (dB)
BITS      = 8            # ADC 位深
MAX_VAL   = 2**BITS - 1  # = 255

# 传感器积分时延（关键修正参数）
n_per     = OVERFS // FS              # = 25（每传感器采样对应的过采样点数）
INT_DELAY = (n_per - 1) / (2 * OVERFS)  # = 2.4ms

# 10个波长通道（nm）
WLS10 = [450, 494, 539, 583, 628, 672, 717, 761, 806, 850]

# 5个LED频率（Hz）× 2正交相位 = 10通道
FREQS  = [13, 23, 33, 43, 53]
PHASES = [0, np.pi/2]   # I通道=0°，Q通道=90°

# 通道映射：(频率, 相位, 波长)
PAIRS = [(FREQS[k//2], PHASES[k%2], WLS10[k]) for k in range(10)]

# ─────────────────────────────────────────────────────
# 2. 传感器与 LED 光谱模型
# ─────────────────────────────────────────────────────

wavelengths = np.arange(400, 1001, 1.0)  # 光谱轴 (nm)

def sensor_response(wl):
    """硅基传感器光谱响应（归一化，400-1000nm）"""
    r = np.zeros_like(wl, dtype=float)
    m = (wl >= 400) & (wl <= 1000)
    wm = wl[m]
    r[m] = np.clip((wm-380)/120, 0, 1) * np.clip((1020-wm)/300, 0, 1)
    return r / r.max()

sr = sensor_response(wavelengths)

def led_spectrum(center_wl, fwhm=20):
    """LED 光谱（高斯线型，半峰全宽 fwhm nm）"""
    return np.exp(-4*np.log(2)*((wavelengths - center_wl)/fwhm)**2)

def compute_weights(wl_list, scale=0.2):
    """LED-传感器有效权重：光谱加权积分，归一化至最大值=scale"""
    raw = np.array([
        np.trapezoid(led_spectrum(wl) * sr, wavelengths)
        for wl in wl_list
    ])
    return raw / raw.max() * scale

ws10 = compute_weights(WLS10)

# ─────────────────────────────────────────────────────
# 3. 目标反射率模型（仿真用时变信号）
# ─────────────────────────────────────────────────────

def target_reflectance(t, wl):
    """
    仿真目标物体在波长 wl 处的时变反射率。
    基础值随波长高斯分布，叠加两个频率分量。
    变化频率：f1 = 0.5+(wl-450)/400 Hz，f2 = 1.7×f1
    """
    f1 = 0.5 + (wl - 450) / 400
    f2 = 1.7 * f1
    r  = (0.4 + 0.3 * np.exp(-((wl-600)/150)**2)
          + 0.15 * np.sin(2*np.pi*f1*t)
          + 0.08 * np.sin(2*np.pi*f2*t + np.pi/3))
    return np.clip(r, 0.05, 0.95)

# ─────────────────────────────────────────────────────
# 4. 传感器信号合成（仿真）
# ─────────────────────────────────────────────────────

def synthesize_sensor_signal(true_refl, seed=42):
    """
    合成传感器信号：
      1. 生成各 LED 的 PWM 调制贡献（高采样率）
      2. 传感器积分采样（平均 n_per 点）
      3. 叠加高斯电子噪声（SNR=42dB）
      4. 8-bit 量化
    """
    np.random.seed(seed)
    t_hi = np.linspace(0, T_SIM, int(T_SIM*OVERFS), endpoint=False)
    ns   = int(T_SIM * FS)

    sig_continuous = np.zeros(len(t_hi))
    for f, phi, wl in PAIRS:
        iw = WLS10.index(wl)
        pwm = (signal.square(2*np.pi*f*t_hi + phi, duty=0.5) + 1) / 2
        sig_continuous += ws10[iw] * pwm * true_refl[wl]

    # 传感器积分
    sig_sampled = sig_continuous[:ns*n_per].reshape(ns, n_per).mean(axis=1)

    # 电子噪声 + 量化
    sig_rms   = np.sqrt(np.mean(sig_sampled**2))
    noise_std = sig_rms / 10**(SNR_DB / 20)
    noise     = np.random.normal(0, noise_std, ns)
    sig_q     = np.clip(
        np.round((sig_sampled + noise) * MAX_VAL),
        0, MAX_VAL
    ) / MAX_VAL

    t_sensor = np.linspace(0, T_SIM, ns, endpoint=False)
    return sig_q, t_sensor, t_hi

# ─────────────────────────────────────────────────────
# 5. IQ 锁相解调（核心算法）
# ─────────────────────────────────────────────────────

def iq_phase_lockin(sig_q, f, ws_val, phi, t_sensor,
                    fs=FS, overfs=OVERFS, cutoff=CUTOFF):
    """
    IQ 同步解调：从叠加信号中提取单个通道的反射率。

    参数
    ----
    sig_q    : 传感器量化信号
    f        : LED 调制频率 (Hz)
    ws_val   : 该通道的传感器权重
    phi      : LED 初始相位（I通道=0，Q通道=π/2）
    t_sensor : 采样时间轴

    原理
    ----
    1. 生成相位修正参考：ref = sin/cos(2πf(t+delay))
       delay 补偿传感器积分引入的时延 τ=(n_per-1)/(2×OVERFS)
    2. 混频 + LPF：lp = LPF[sig × ref]
    3. 校准：R = π×lp / ws_val
       （由 PWM d=0.5 基频分量 = 2/π × sin/cos，混频后 LPF = ws×R/π）
    """
    n_per_ = overfs // fs
    delay  = (n_per_ - 1) / (2 * overfs)

    if phi < 0.1:  # I通道
        ref = np.sin(2 * np.pi * f * (t_sensor + delay))
    else:           # Q通道
        ref = np.cos(2 * np.pi * f * (t_sensor + delay))

    sos = butter(4, cutoff/(fs/2), btype='low', output='sos')
    lp  = sosfiltfilt(sos, sig_q * ref)
    return np.pi * lp / ws_val

def reconstruct_all(sig_q, t_sensor, trim=1.0):
    """对所有10个通道执行 IQ 解调，返回有效区间的重建反射率"""
    vl = (t_sensor >= trim) & (t_sensor <= T_SIM - trim)
    recon = {}
    for f, phi, wl in PAIRS:
        iw = WLS10.index(wl)
        R = iq_phase_lockin(sig_q, f, ws10[iw], phi, t_sensor)
        recon[wl] = R[vl]
    t_valid = t_sensor[vl]
    return recon, t_valid

# ─────────────────────────────────────────────────────
# 6. 性能评估
# ─────────────────────────────────────────────────────

def evaluate(recon, t_valid, true_refl, t_hi):
    """计算各通道的 RMSE 和皮尔逊相关系数"""
    results = {}
    for wl in WLS10:
        ri   = np.interp(t_valid, t_hi, true_refl[wl])
        r    = recon[wl]
        rmse = np.sqrt(np.mean((r - ri)**2))
        corr = np.corrcoef(r, ri)[0, 1]
        results[wl] = {'rmse': rmse, 'corr': corr}
    return results

# ─────────────────────────────────────────────────────
# 7. 主程序
# ─────────────────────────────────────────────────────

if __name__ == '__main__':
    np.random.seed(42)

    # 生成真实反射率
    t_hi = np.linspace(0, T_SIM, int(T_SIM*OVERFS), endpoint=False)
    true_refl = {wl: target_reflectance(t_hi, wl) for wl in WLS10}

    # 合成传感器信号
    print('合成传感器信号...')
    sig_q, t_sensor, t_hi_out = synthesize_sensor_signal(true_refl)

    # 重建
    print('IQ 锁相解调重建...')
    recon, t_valid = reconstruct_all(sig_q, t_sensor)

    # 评估
    metrics = evaluate(recon, t_valid, true_refl, t_hi)

    print(f'\n{"通道":<10} {"LED频率":>8} {"相位":>6} {"RMSE":>10} {"皮尔逊r":>10}')
    print('-' * 50)
    for k, (f, phi, wl) in enumerate(PAIRS):
        ph_deg = 0 if phi < 0.1 else 90
        m = metrics[wl]
        print(f'{wl}nm     {f:>6}Hz {"I" if ph_deg==0 else "Q":>6}'
              f'  {m["rmse"]:>10.5f}  {m["corr"]:>10.5f}')
    mean_rmse = np.mean([metrics[wl]['rmse'] for wl in WLS10])
    print(f'{"均值":<10} {"":>8} {"":>6}  {mean_rmse:>10.5f}')

    # 基本可视化
    fig, axes = plt.subplots(5, 2, figsize=(14, 18), sharex=True)
    for k, (f, phi, wl) in enumerate(PAIRS):
        ax = axes[k//2, k%2]
        ri = np.interp(t_valid, t_hi, true_refl[wl])
        ax.plot(t_valid, ri, lw=2, alpha=0.6, label='真实值')
        ax.plot(t_valid, recon[wl], 'k--', lw=1.5,
                label=f'重建 RMSE={metrics[wl]["rmse"]:.4f}')
        ax.set_ylim(0, 1.1)
        ax.set_title(f'{wl}nm | {f}Hz-{"I" if phi<0.1 else "Q"}', fontsize=10)
        ax.legend(fontsize=8)
    axes[-1, 0].set_xlabel('时间(s)'); axes[-1, 1].set_xlabel('时间(s)')
    fig.suptitle('IQ相位编码系统 — 10通道重建结果', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig('reconstruction_result.png', dpi=150, bbox_inches='tight')
    plt.show()
    print('\n结果图已保存: reconstruction_result.png')
