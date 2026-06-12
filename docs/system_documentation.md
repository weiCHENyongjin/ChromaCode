# 多波长LED相位编码光强采集与反射率重建系统

**版本**: v1（第一阶段）  
**日期**: 2026-06-05  
**传感器**: 200 Hz / 8-bit / SNR≈42 dB  
**通道数**: 10（5频率 × 2正交相位）

---

## 1. 数学原理

### 1.1 PWM 方波的傅里叶分解

对于频率 $f$、占空比 $d$、初始相位 $\varphi$ 的 PWM 方波信号（高电平为 1，低电平为 0），其傅里叶级数为：

$$
P(t; f, d, \varphi) = d + \sum_{n=1}^{\infty} \frac{2\sin(n\pi d)}{n\pi} \cos\!\bigl(2\pi n f t - n\pi d + n\varphi\bigr)
$$

**关键性质**：

- 当 $d = 0.5$（50% 占空比）时，**所有偶次谐波系数为零**（$\sin(n\pi \times 0.5) = 0$ 对偶数 $n$），仅保留奇次谐波。  
- 基频（$n=1$）幅值为 $\dfrac{2\sin(\pi d)}{\pi}$，在 $d=0.5$ 时取得最大值 $\dfrac{2}{\pi} \approx 0.637$。

### 1.2 正交相位编码（IQ 调制）

令同一频率 $f$ 的两路 LED 的初始相位分别为 $\varphi_1 = 0$ 和 $\varphi_2 = \pi/2$：

$$
P_I(t) = P(t;\, f,\, 0.5,\, 0) = 0.5 + \frac{2}{\pi}\sin(2\pi ft) + \cdots
$$

$$
P_Q(t) = P(t;\, f,\, 0.5,\, \tfrac{\pi}{2}) = 0.5 + \frac{2}{\pi}\cos(2\pi ft) + \cdots
$$

两路 LED 分别照射不同波长的目标，目标对两波长的反射率分别为 $R_I(t)$ 和 $R_Q(t)$。传感器（光谱加权积分权重分别为 $w_I,\, w_Q$）采集到的叠加信号为：

$$
s(t) = w_I R_I(t) P_I(t) + w_Q R_Q(t) P_Q(t) + \text{（其他频率组贡献）} + \epsilon(t)
$$

### 1.3 锁相检测（Lock-in Detection）

以频率 $f$ 为参考，对传感器信号做 **同步解调**：

$$
\tilde{I}(t) = \mathrm{LPF}\!\bigl[s(t) \cdot \sin(2\pi f t)\bigr], \quad
\tilde{Q}(t) = \mathrm{LPF}\!\bigl[s(t) \cdot \cos(2\pi f t)\bigr]
$$

其中 $\mathrm{LPF}$ 为低通滤波器（截止频率 $f_c$）。  
利用正交关系 $\mathrm{LPF}[\sin^2] = 1/2$、$\mathrm{LPF}[\cos^2] = 1/2$、$\mathrm{LPF}[\sin\cos] = 0$，代入后得：

$$
\tilde{I}(t) \approx \frac{w_I R_I(t)}{\pi}, \qquad \tilde{Q}(t) \approx \frac{w_Q R_Q(t)}{\pi}
$$

**反射率重建公式**：

$$
\boxed{R_I(t) = \frac{\pi \,\tilde{I}(t)}{w_I}, \qquad R_Q(t) = \frac{\pi \,\tilde{Q}(t)}{w_Q}}
$$

### 1.4 传感器积分时延的相位修正

实际传感器对每个采样点进行 **有限时间积分**（曝光），采样频率为 $F_s$、过采样仿真频率为 $F_\text{sim}$：

$$
T_\text{int} = \frac{1}{F_s}, \qquad
n_\text{per} = \frac{F_\text{sim}}{F_s}
$$

积分操作等效于对连续信号施加 boxcar 均值滤波，其中心时延为：

$$
\tau = \frac{n_\text{per} - 1}{2 \cdot F_\text{sim}}
$$

该时延使 LED 的载波相对于参考信号产生相移 $\Delta\phi_f = 2\pi f \tau$，导致高频 LED 的锁相输出幅值严重衰减（如 97 Hz 时 $\cos(\Delta\phi) \approx 0.09$，即只剩 9% 幅值）。

**修正方法**：在参考信号中加入时延补偿：

$$
\text{ref}_I(t) = \sin\!\bigl(2\pi f(t + \tau)\bigr), \quad
\text{ref}_Q(t) = \cos\!\bigl(2\pi f(t + \tau)\bigr)
$$

### 1.5 时间分辨率

锁相检测的时间分辨率由低通滤波器截止频率 $f_c$ 决定，可追踪的目标信号最高频率为 $f_c$（4阶 Butterworth 在 $f_c$ 处约 −3 dB）。在本方案中：

$$
f_c = 3.5 \text{ Hz} \implies \text{可追踪信号} \leq 3.5 \text{ Hz}
$$

### 1.6 LED 频率设计约束

为使不同频率组之间的串扰最小，需满足：

1. **基频拍频隔离**：$|f_i - f_j| > 2f_c = 7 \text{ Hz}$（所有频率对）  
2. **谐波拍频隔离**：$\bigl|\mathrm{alias}(k f_i, F_s) - f_j\bigr| > 2f_c$（$k = 3, 5, 7, \ldots$）  
3. **奈奎斯特限制**：$f_{\max} < F_s/2 = 100 \text{ Hz}$

**最优间距的折中**：实验扫描（见第4节）表明，频率间距 $G \approx 10 \text{ Hz}$ 时综合 RMSE 最小，其最高 LED 频率为 53 Hz，使所有谐波混叠效应和相位修正误差均处于可接受范围内。

---

## 2. 算法说明

### 2.1 系统架构

```
目标物体（时变反射率）
      │
  LED照明（5频率×2相位=10路）
      │
  单像素传感器（积分采样）
      │
  ADC量化（8-bit）
      │
  IQ锁相解调（5次×2路）
      │
  低通滤波（4阶Butterworth, 3.5Hz）
      │
  幅值校准（÷ ws_k/π）
      │
  10路时变反射率序列
```

### 2.2 传感器权重计算

传感器对波长 $\lambda$ 的响应 $\eta(\lambda)$（近似硅基传感器曲线），LED $k$ 的有效权重为光谱加权积分：

$$
w_k = \int_{400}^{1000} I_k(\lambda)\, \eta(\lambda)\, d\lambda
$$

其中 $I_k(\lambda) = \exp\!\left[-\frac{4\ln 2 \cdot (\lambda - \lambda_k^0)^2}{\Delta\lambda_k^2}\right]$（高斯线型，半峰全宽 $\Delta\lambda_k = 20 \text{ nm}$）。

所有权重归一化至最大值 0.2（使最大信号幅值约占满量程的 20%）。

### 2.3 仿真信号合成

1. 以 $F_\text{sim} = 5000 \text{ Hz}$ 过采样生成连续时间信号
2. 对每个 LED（频率 $f_k$，相位 $\varphi_k \in \{0, \pi/2\}$，占空比 0.5）：

$$
s_k(t) = w_k \cdot \frac{P_k(t) + 1}{2} \cdot R_k(t)
$$

3. 叠加全部 LED 贡献：$s(t) = \sum_k s_k(t)$
4. 传感器积分（对 $n_\text{per} = F_\text{sim}/F_s = 25$ 个过采样点取均值）
5. 加入高斯电子噪声（$\sigma_e = s_\text{rms} / 10^{\text{SNR}/20}$）
6. 8-bit 量化：$s_q[n] = \mathrm{round}(s[n] \times 255) / 255$

### 2.4 IQ 锁相重建

核心解调（每通道一条按相位匹配的参考信号）的概念实现如下；完整实现见
`python/spectral_reconstruction.py` 中的 `SpectralReconstructor.reconstruct`（I 通道
$\varphi=0$ 时 $\sin$ 即标准正弦，Q 通道 $\varphi=\pi/2$ 时 $\sin(\cdot+\pi/2)=\cos$，
二者由同一公式统一表达）：

```python
# 每个通道：频率 f、相位 phi、权重 w
delay = integration_time / 2          # 传感器积分群延迟（见 §1.4 / 理论文档 §6）
ref   = np.sin(2 * np.pi * f * (t_sensor + delay) + phi)   # 相位补偿的匹配参考

sos   = butter(4, cutoff / (fs / 2), btype='low', output='sos')
lp    = sosfiltfilt(sos, sig_q * ref)  # 4 阶 Butterworth，零相位双向滤波

R     = np.pi * lp / w                 # 幅值校准 → 反射率
```

**校准公式推导**：对于 $d=0.5$ 的 PWM，混频后低通分量为：

$$
\tilde{I} = \mathrm{LPF}[s \cdot \sin_\text{ref}] = \frac{w_I R_I}{\pi}
\implies R_I = \frac{\pi \tilde{I}}{w_I}
$$

---

## 3. 实验配置与结果

### 3.1 参数配置

| 参数 | 值 |
|------|---|
| 传感器采样率 $F_s$ | 200 Hz |
| 过采样仿真率 $F_\text{sim}$ | 5000 Hz |
| ADC 位深 | 8-bit（255 级）|
| 信噪比 SNR | 42 dB |
| 锁相低通截止 $f_c$ | 3.5 Hz |
| 低通滤波器阶数 | 4 阶 Butterworth（双向，零相位）|
| 仿真时长 | 12 s |
| 有效评估区间 | 1 s ～ 11 s（去掉首尾边缘效应）|

### 3.2 LED 配置（10通道最优方案）

频率间距 $G = 10 \text{ Hz}$，频率组 $[13, 23, 33, 43, 53] \text{ Hz}$。

| 频率组 | I 通道（phase=0°） | Q 通道（phase=90°）|
|--------|-------------------|-------------------|
| 13 Hz  | 450 nm | 494 nm |
| 23 Hz  | 539 nm | 583 nm |
| 33 Hz  | 628 nm | 672 nm |
| 43 Hz  | 717 nm | 761 nm |
| 53 Hz  | 806 nm | 850 nm |

**LED 光谱参数**：中心波长如上，FWHM = 20 nm，高斯线型，强度归一化，占空比 = 50%。

### 3.3 目标反射率模型（仿真用）

每个波长 $\lambda$ 的时变反射率：

$$
R(t;\lambda) = \text{clip}\!\left[\,0.4 + 0.3 e^{-\left(\frac{\lambda-600}{150}\right)^2}
+ 0.15\sin(2\pi f_1 t) + 0.08\sin(2\pi f_2 t + \tfrac{\pi}{3}),\; 0.05,\; 0.95\right]
$$

其中 $f_1 = 0.5 + \dfrac{\lambda - 450}{400} \text{ Hz}$，$f_2 = 1.7 f_1$。

各通道信号变化频率范围：$f_2 \in [0.85, 2.55] \text{ Hz}$（最快变化为 850 nm 通道）。

### 3.4 重建结果

| 通道 | LED 频率 | 相位 | 信号变化频率 | RMSE | 皮尔逊 $r$ |
|------|---------|------|------------|------|-----------|
| 450 nm | 13 Hz | I | 0.50 + 0.85 Hz | **0.0246** | **0.979** |
| 494 nm | 13 Hz | Q | 0.61 + 1.04 Hz | **0.0168** | **0.991** |
| 539 nm | 23 Hz | I | 0.72 + 1.23 Hz | **0.0212** | **0.993** |
| 583 nm | 23 Hz | Q | 0.83 + 1.42 Hz | **0.0232** | **0.995** |
| 628 nm | 33 Hz | I | 0.95 + 1.61 Hz | 0.0556 | 0.930 |
| 672 nm | 33 Hz | Q | 1.06 + 1.79 Hz | 0.0466 | 0.961 |
| 717 nm | 43 Hz | I | 1.17 + 1.98 Hz | 0.0633 | 0.931 |
| 761 nm | 43 Hz | Q | 1.28 + 2.17 Hz | 0.0672 | 0.895 |
| 806 nm | 53 Hz | I | 1.39 + 2.36 Hz | **0.0493** | **0.986** |
| 850 nm | 53 Hz | Q | 1.50 + 2.55 Hz | 0.0717 | 0.981 |
| **均值** | | | | **0.044** | **0.964** |

---

## 4. 结果分析

### 4.1 相位编码的串扰隔离

对于 $d = 0.5$ 的 PWM，奇次谐波的相位随 LED 初始相位线性偏移（$n$ 次谐波的相位偏移为 $n\varphi$）。因此，相位为 $0$ 和 $\pi/2$ 的两路 LED：

- **基频正交**：$\sin(2\pi ft)$ 与 $\cos(2\pi ft)$ 完全正交，I/Q 通道零串扰
- **奇次谐波**：相位偏移为 $n\pi/2$，与参考信号的乘积 LPF 后亦趋近于零
- **偶次谐波全消**：$d=0.5$ 时 $\sin(n\pi \times 0.5) = 0$（$n$ 为偶数），彻底消除偶次谐波干扰

这使本方案远优于使用 $d = 0.25/0.75$ 的占空比编码方案（后者存在较大偶次谐波）。

### 4.2 传感器积分时延的影响

| LED 频率 | 相位误差（未修正）| 幅值保留率（未修正）| 修正后效果 |
|---------|---------------|------------------|---------|
| 13 Hz | 11.2° | 98.1% | 无需关注 |
| 23 Hz | 19.8° | 94.1% | 已修正 |
| 33 Hz | 28.4° | 87.9% | 已修正 |
| 43 Hz | 37.0° | 79.9% | 已修正 |
| 53 Hz | 45.6° | 70.2% | 已修正 |

时延 $\tau = 12/5000 = 2.4 \text{ ms}$，相移 $\Delta\phi = 2\pi f \cdot 0.0024 \text{ rad}$。

### 4.3 频率间距优化

对频率间距 $G \in \{6, 8, 10, 12, 14, 16, 18, 20, 22\}$ Hz 进行实验扫描：

| 间距 $G$ | 频率组 | 均值 RMSE | 特征 |
|---------|-------|---------|------|
| 6 Hz | [13,19,25,31,37] | 0.066 | 短程，高频串扰 |
| 8 Hz | [13,21,29,37,45] | 0.053 | 改善 |
| **10 Hz** | **[13,23,33,43,53]** | **0.044** | **最优** |
| 12 Hz | [13,25,37,49,61] | 0.064 | 628nm 谐波冲突 |
| 16 Hz | [13,29,45,61,77] | 0.058 | 较好 |
| 22 Hz | [5,27,49,71,93] | 0.175 | 5Hz×5=25Hz 精确冲突 |

**最优 $G = 10 \text{ Hz}$ 的物理解释**：

1. 最高 LED 频率 53 Hz < 100 Hz（奈奎斯特），且相位误差约 45°（可以被修正）
2. 相邻频率间距 10 Hz > $2f_c = 7 \text{ Hz}$，满足基频隔离
3. 3次谐波最小距离 6 Hz（如 $13\times3=39$ Hz，距 33 Hz 为 6 Hz），处于可接受水平

### 4.4 与原始5通道方案的对比

| 指标 | 原始方案 | 本方案（IQ编码）|
|------|---------|---------------|
| 通道数 | 5 | **10** |
| 传感器采样率 | 500 Hz | **200 Hz** |
| ADC 位深 | 12-bit | **8-bit** |
| LED 数量 | 5 | **10**（但占2倍带宽×2倍相位）|
| 均值 RMSE（5通道对比） | 0.056 | **0.044**（10通道）|
| 可追踪信号频率 | ≤0.5 Hz（STFT限制）| **≤3.5 Hz**（LPF限制）|

### 4.5 局限性分析

1. **快变通道精度下降**：717–850 nm 通道（$f_2 = 2.0$–$2.55$ Hz）的 RMSE（0.063–0.072）高于慢变通道，因为信号频率接近 LPF 截止频率 3.5 Hz（幅值衰减约 −3 dB）。

2. **谐波混叠残余**：$53 \times 3 = 159 \text{ Hz} \to 41 \text{ Hz}$（混叠），距 43 Hz 仅 2 Hz，使 628/717 nm 通道存在轻微系统误差（RMSE≈0.05–0.06）。

3. **IIR 滤波器边缘效应**：低截止频率（3.5 Hz）的双向 Butterworth 滤波器在信号边缘产生约 1 s 的过渡段，需去掉首尾各 1 s。

---

## 5. 总结

### 5.1 核心贡献

本文提出并验证了一种基于**正交相位编码**（IQ Phase Encoding）的多波长 LED 傅里叶分析系统：

- **每个 LED 频率携带 2 个正交通道**，5 个频率组实现 10 通道，通道密度翻倍
- **相位编码（$\varphi = 0$ / $\pi/2$）而非占空比编码**，因 $d = 0.5$ 时偶次谐波完全消除，串扰更低
- **锁相检测而非 STFT**，时间分辨率由 LPF 截止频率（3.5 Hz）决定，可追踪更快信号
- **传感器积分时延相位修正**，解决了高频 LED（53 Hz）在 200 Hz 采样率下的相位失真问题

### 5.2 最优系统参数

```
传感器:  Fs = 200 Hz,  8-bit ADC,  SNR ≥ 42 dB
LED:     [13, 23, 33, 43, 53] Hz（间距 G = 10 Hz）
         各对相位: 0° (I通道) 和 90° (Q通道)
         占空比: 50%（消除偶次谐波）
重建:    锁相解调 + 4阶Butterworth LPF（截止3.5 Hz）
         相位修正延迟 τ = 12/5000 = 2.4 ms
性能:    10通道均值 RMSE = 0.044，皮尔逊 r = 0.964
```

### 5.3 后续研究方向

1. **实物验证**：将仿真中的 PWM 相位控制、传感器积分时延测量移植到硬件
2. **更多通道**：在 200 Hz 以内加入更多频率组，或通过占空比+相位组合编码达到 15+ 通道
3. **自适应校准**：通过白参考板自动标定各通道的 $w_k$ 权重，消除光源老化影响
4. **提高 Fs**：采样率从 200 Hz 提高至 1000 Hz 可彻底消除谐波混叠，850 nm 通道精度预计提升至 RMSE ≈ 0.01

---

## 附录：源代码

为保持"单一事实来源"（避免文档内嵌代码与实现脱节），完整可复现代码不再内联于此，
请直接查阅仓库源码（均含详细注释，遵循 `../CODE_STYLE.md` 的统一规范）：

| 文件 | 内容 |
|------|------|
| [`../python/iq_sensing_system.py`](../python/iq_sensing_system.py) | 前向仿真器（`ForwardSimulator`）、目标反射率模型、主程序 |
| [`../python/spectral_reconstruction.py`](../python/spectral_reconstruction.py) | 锁相重建核心（`SpectralReconstructor`）、光谱模型、配置加载 |
| [`../python/make_figures.py`](../python/make_figures.py) | 中英双语结果图生成 |
| [`../config/example_config.yaml`](../config/example_config.yaml) | 传感器 + LED 配置模板 |
| [`../python/examples/example_usage.py`](../python/examples/example_usage.py) | 最小可运行示例 |
| [`../python/tests/test_reconstruction.py`](../python/tests/test_reconstruction.py) | 端到端等价性测试（复现均值 RMSE ≈ 0.044） |

运行：

```bash
pip install -r ../python/requirements.txt
python ../python/iq_sensing_system.py        # 仿真 + 重建 + 结果图
python ../python/tests/test_reconstruction.py     # 验证精度
```

完整数学理论（定义、引理、定理及证明）见 [`mathematical_theory.md`](mathematical_theory.md)。

---

*文档结束*
