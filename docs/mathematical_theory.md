# 多波长 LED 相位编码光强采集与多光谱重建——数学理论

**Mathematical Theory of IQ Phase-Encoded Multi-Wavelength Reflectance Reconstruction**

> 本文给出该系统从基础假设到重建公式的完整数学理论，包括定义、引理、定理及其证明、误差分析与设计约束的推导。
> 记号与符号见 §0。配套实现见 [`../code/spectral_reconstruction.py`](../code/spectral_reconstruction.py)，工程化文档见 [`system_documentation.md`](system_documentation.md)。
>
> *GitHub 已支持 LaTeX 数学渲染；若本地阅读器不渲染公式，请用支持 KaTeX/MathJax 的查看器。*

---

## 目录

- §0 记号与基本假设
- §1 基础理论：PWM 方波的傅里叶分解
- §2 正交关系引理
- §3 前向测量模型
- §4 锁相解调（单通道重建定理）
- §5 IQ 正交相位编码（通道分离定理与串扰分析）
- §6 传感器积分时延与相位补偿
- §7 多频串扰与频率规划约束
- §8 时间分辨率
- §9 光谱加权、波长重采样与参考标定
- §10 误差分析（噪声传播与量化）
- §11 定理汇总与设计准则
- 参考文献

---

## §0 记号与基本假设

### 0.1 符号表

| 符号 | 含义 |
|------|------|
| $K$ | 通道总数（本方案 $K=10$） |
| $\lambda_k$ | 第 $k$ 通道的中心波长 (nm) |
| $f_m$ | 第 $m$ 个 LED 调制频率 (Hz)，$m=1,\dots,M$（本方案 $M=5$） |
| $d$ | PWM 占空比（duty cycle），本方案 $d=\tfrac12$ |
| $\varphi_k$ | 第 $k$ 通道 LED 的初始相位（I 通道 $0$，Q 通道 $\pi/2$） |
| $R_k(t)$ | 目标在 $\lambda_k$ 处的时变反射率，$R_k:\mathbb{R}\to[0,1]$ |
| $w_k$ | 第 $k$ 通道的 LED–传感器有效权重（光谱重叠积分） |
| $F_s$ | 传感器采样率 (Hz)，本方案 $200$ |
| $T_\mathrm{int}$ | 传感器积分（曝光）时间 (s) |
| $f_c$ | 锁相低通滤波器截止频率 (Hz)，本方案 $3.5$ |
| $\varepsilon(t)$ | 加性电子噪声 |
| $\eta(\lambda)$ | 传感器光谱响应（可选，常未知） |
| $I_k(\lambda)$ | 第 $k$ 通道 LED 的发射光谱 |
| $\langle\cdot\rangle$ | 长时间平均算子 $\langle g\rangle=\lim_{T\to\infty}\tfrac1T\int_0^T g$ |
| $\mathrm{LPF}[\cdot]$ | 截止 $f_c$ 的零相位低通滤波（局部平均的物理实现） |

### 0.2 基本假设

- **(A1) 慢变假设**：每个反射率 $R_k(t)$ 的带宽 $B_k$ 满足 $B_k \le f_c \ll f_m$。即反射率在一个调制周期 $1/f_m$ 内近似为常数。
- **(A2) 线性叠加**：传感器输出与各 LED 贡献线性叠加，光电响应在工作区间线性。
- **(A3) 加性噪声**：电子噪声 $\varepsilon(t)$ 为零均值、与信号独立的广义平稳过程。
- **(A4) 频率规划**：调制频率 $\{f_m\}$ 满足 §7 的隔离条件，使跨频串扰落在 $\mathrm{LPF}$ 阻带内。

---

## §1 基础理论：PWM 方波的傅里叶分解

### 定义 1.1（单极性 PWM 脉冲）

给定频率 $f$、占空比 $d\in(0,1)$、初始相位 $\varphi$，定义相位变量 $\theta(t)=2\pi f t+\varphi$。单极性 PWM 脉冲为

$$
p(t;f,d,\varphi)=\begin{cases}1, & \theta(t)\bmod 2\pi \in [0,\,2\pi d),\\[2pt]0, & \text{否则.}\end{cases}
$$

它是周期 $T=1/f$ 的方波，取值 $\{0,1\}$。（这正是代码中 `(scipy.signal.square(θ, duty=d)+1)/2`。）

### 定理 1.2（PWM 傅里叶级数）

$$
\boxed{\;p(t;f,d,\varphi)=d+\sum_{n=1}^{\infty}\frac{2\sin(n\pi d)}{n\pi}\,\cos\!\big(2\pi n f t+n\varphi-n\pi d\big).\;}
$$

**证明.** 以 $\theta$ 为自变量，$p$ 是 $2\pi$ 周期函数，复傅里叶系数

$$
c_n=\frac{1}{2\pi}\int_{0}^{2\pi d} e^{-in\theta}\,d\theta
=\frac{1}{2\pi}\cdot\frac{1-e^{-i\,2\pi n d}}{in},\qquad n\neq 0,
$$

且 $c_0=\frac{1}{2\pi}\int_0^{2\pi d}d\theta=d$。利用

$$
1-e^{-i2\pi nd}=e^{-i\pi nd}\big(e^{i\pi nd}-e^{-i\pi nd}\big)=e^{-i\pi nd}\cdot 2i\sin(n\pi d),
$$

得

$$
c_n=\frac{1}{2\pi in}\cdot 2i\sin(n\pi d)\,e^{-i\pi nd}=\frac{\sin(n\pi d)}{n\pi}\,e^{-i\pi nd}.
$$

代回 $p=\sum_n c_n e^{in\theta}=c_0+\sum_{n\ge1}2\,\mathrm{Re}\!\big[c_n e^{in\theta}\big]$，并将 $\theta=2\pi ft+\varphi$ 代入，即得定理所述级数。$\blacksquare$

### 引理 1.3（$d=\tfrac12$：偶次谐波全消，基频取极大）

当 $d=\tfrac12$ 时：

1. 对所有**偶数** $n$，谐波系数为零（$\sin(n\pi/2)=0$）；
2. 基频（$n=1$）化简为

$$
\text{fundamental}=\frac{2}{\pi}\,\sin\!\big(2\pi f t+\varphi\big),
$$

其幅值 $\tfrac{2}{\pi}=\max_d \tfrac{2\sin(\pi d)}{\pi}$ 在 $d=\tfrac12$ 处取得。

**证明.** (1) $\sin(n\pi/2)=0\Leftrightarrow n$ 为偶数。(2) 取 $n=1,d=\tfrac12$：定理 1.2 给出 $\tfrac{2\sin(\pi/2)}{\pi}\cos(2\pi ft+\varphi-\pi/2)=\tfrac{2}{\pi}\cos(2\pi ft+\varphi-\tfrac\pi2)=\tfrac{2}{\pi}\sin(2\pi ft+\varphi)$。基频幅值 $\tfrac{2\sin(\pi d)}{\pi}$ 在 $d=\tfrac12$ 处对 $d$ 求导为零且为极大。$\blacksquare$

> **意义.** 偶次谐波全消是选择 $d=\tfrac12$ 的根本原因：它消除了 $2f,4f,\dots$ 这些最易与相邻频率组冲突的谐波，并使基频能量最大、信噪比最优。相位编码（§5）正是建立在"基频为纯 $\tfrac{2}{\pi}\sin(2\pi ft+\varphi)$"这一干净形式之上。

由引理 1.3，本方案中每个 LED 的有效照明可写为

$$
p_k(t)=\underbrace{\tfrac12}_{\text{直流}}+\underbrace{\tfrac{2}{\pi}\sin(2\pi f_{m(k)}t+\varphi_k)}_{\text{基频}}+\underbrace{\sum_{n\ \text{odd}\ge3}\tfrac{2\sin(n\pi/2)}{n\pi}\cos(\cdots)}_{\text{奇次谐波}}.\tag{1.1}
$$

---

## §2 正交关系引理

### 引理 2.1（正弦正交性）

设 $f,g>0$。在长时间平均 $\langle\cdot\rangle$ 下，

$$
\langle \sin(2\pi f t+\alpha)\,\sin(2\pi g t+\beta)\rangle=
\begin{cases}\tfrac12\cos(\alpha-\beta), & f=g,\\[2pt]0,& f\neq g.\end{cases}
$$

**证明.** 积化和差：$\sin A\sin B=\tfrac12[\cos(A-B)-\cos(A+B)]$。其中 $A-B=2\pi(f-g)t+(\alpha-\beta)$，$A+B=2\pi(f+g)t+(\alpha+\beta)$。当 $f\neq g$，两项均为非零频率余弦，时间平均为 $0$。当 $f=g$，$\cos(A-B)=\cos(\alpha-\beta)$ 为常数，$\cos(A+B)$ 平均为 $0$，故得 $\tfrac12\cos(\alpha-\beta)$。$\blacksquare$

### 推论 2.2（同频正交相位）

取 $f=g$ 且 $\alpha-\beta=\pm\pi/2$，则 $\langle\sin\cdot\sin\rangle=\tfrac12\cos(\pm\pi/2)=0$。这是 I/Q 两路（相差 $90^\circ$）零串扰的代数核心。

### 引理 2.3（慢变信号可移出平均）

在假设 (A1) 下，若 $R(t)$ 带宽 $\le f_c$，而 $u(t)$ 的能量集中在 $\ge f_m-f_c>2f_c$ 的频带，则

$$
\mathrm{LPF}\big[R(t)\,u(t)\big]\approx R(t)\,\mathrm{LPF}[u(t)],
$$

误差为 $O(B_R/f_m)$。

**证明（要点）.** 乘积在频域是卷积 $\hat R * \hat u$。$\hat R$ 支撑于 $[-f_c,f_c]$，$\hat u$ 的谱线在 $\pm f_m$ 附近。卷积仅在 $\hat u$ 自身落入 $[-f_c,f_c]$ 处对 $\mathrm{LPF}$ 输出有贡献，此时 $\hat R$ 近似为其在该窄带上的取值，等价于把 $R(t)$ 视为缓慢包络移出滤波器。带宽比 $B_R/f_m$ 控制残差。$\blacksquare$

---

## §3 前向测量模型

### 3.1 连续叠加信号

第 $k$ 通道 LED 照射波长 $\lambda_k$，目标反射率 $R_k(t)$，传感器有效权重 $w_k$（见 §9）。由 (A2) 叠加，传感器入射光强（连续时间）为

$$
x(t)=\sum_{k=1}^{K} w_k\,R_k(t)\,p_k(t).\tag{3.1}
$$

### 3.2 积分采样算子

实际传感器对每个采样点在曝光窗内积分。定义**前向 boxcar 积分采样算子** $\mathcal{S}$：

$$
(\mathcal{S}x)[n]=\frac{1}{T_\mathrm{int}}\int_{nT_s}^{nT_s+T_\mathrm{int}}x(\tau)\,d\tau,\qquad T_s=\frac{1}{F_s},\ \ n=0,1,\dots\tag{3.2}
$$

其离散实现（仿真中 $F_\mathrm{sim}=$`OVERFS`，每采样 $n_\mathrm{per}=F_\mathrm{sim}/F_s$ 个过采样点取均值）为

$$
(\mathcal{S}x)[n]=\frac{1}{n_\mathrm{per}}\sum_{j=0}^{n_\mathrm{per}-1}x\!\big(nT_s+jT_\mathrm{sim}\big),\qquad T_\mathrm{sim}=1/F_\mathrm{sim}.
$$

该算子的相位/幅值效应见 §6。

### 3.3 噪声与量化

加性电子噪声后经 $B$ 位均匀量化（步长 $q=2^{-B}$，满量程归一）：

$$
s[n]=\mathcal{Q}_B\!\Big((\mathcal{S}x)[n]+\varepsilon[n]\Big),\qquad
\mathcal{Q}_B(u)=\frac{\mathrm{round}(u\cdot(2^B-1))}{2^B-1}.\tag{3.3}
$$

噪声标准差由信噪比设定：$\sigma_\varepsilon=\mathrm{RMS}(\mathcal{S}x)\cdot 10^{-\mathrm{SNR}/20}$。量化等效为附加噪声（§10）。$s[n]$ 即"采集到的时序信号"，是重建算法的唯一输入。

---

## §4 锁相解调：单通道重建定理

暂忽略积分时延（§6 补偿后等价）。对目标通道 $m$，构造**匹配参考**

$$
r_m(t)=\sin(2\pi f_m t+\varphi_m).\tag{4.1}
$$

定义解调输出 $y_m(t)=\mathrm{LPF}\big[s(t)\,r_m(t)\big]$。

### 定理 4.1（单通道反射率重建）

在假设 (A1)–(A4) 下，

$$
y_m(t)=\frac{w_m}{\pi}R_m(t)+\text{（阻带残差）},
\qquad\Longrightarrow\qquad
\boxed{\;\hat R_m(t)=\frac{\pi\,y_m(t)}{w_m}=\frac{\pi}{w_m}\,\mathrm{LPF}\big[s\,r_m\big].\;}
$$

**证明.** 由 (3.1)、(1.1)，$s(t)=\sum_k w_k R_k(t)p_k(t)+\varepsilon$。逐项乘 $r_m$ 并取 $\mathrm{LPF}$（即局部 $\langle\cdot\rangle$）：

1. **目标基频项**：$w_m R_m(t)\cdot\tfrac{2}{\pi}\sin(2\pi f_m t+\varphi_m)\cdot\sin(2\pi f_m t+\varphi_m)$。由引理 2.3 移出 $R_m$，再用引理 2.1（$f=g,\alpha=\beta$）得 $\langle\sin^2\rangle=\tfrac12$，贡献 $w_m R_m\cdot\tfrac{2}{\pi}\cdot\tfrac12=\dfrac{w_m R_m}{\pi}$。
2. **目标直流项** $\tfrac12 w_m R_m\cdot r_m$：频率为 $f_m$，落在 $\mathrm{LPF}$ 阻带（$f_m\gg f_c$），平均为 $0$。
3. **目标奇次谐波**（$3f_m,5f_m,\dots$）乘 $r_m$（$f_m$）：频率 $\neq f_m$，引理 2.1 给 $0$。
4. **他通道** $k\neq m$：若 $f_{m(k)}\neq f_m$，基频/谐波与 $r_m$ 频率不等，引理 2.1 给 $0$（跨频残差由 (A4) 保证落入阻带，详见 §7）；若 $f_{m(k)}=f_m$（同频 IQ 配对），由推论 2.2，$90^\circ$ 相位差使其平均为 $0$（详见 §5）。
5. **噪声** $\varepsilon$：零均值，乘有界参考后经窄带 $\mathrm{LPF}$，期望为 $0$（方差见 §10）。

合计 $y_m=\dfrac{w_m}{\pi}R_m+$（阻带与噪声残差）。解出 $\hat R_m$ 即证。$\blacksquare$

> **校准常数 $\pi$ 的来源**：$d=\tfrac12$ 基频幅值 $\tfrac{2}{\pi}$ 与 $\langle\sin^2\rangle=\tfrac12$ 相乘得 $\tfrac1\pi$，故乘 $\pi$ 还原。代码 `np.pi * lp / ws_val` 即此式。

---

## §5 IQ 正交相位编码：通道分离定理与串扰分析

设同一频率 $f$ 承载两路：I 通道（$\varphi_I=0$，反射率 $R_I$，权重 $w_I$）与 Q 通道（$\varphi_Q=\pi/2$，$R_Q,w_Q$）。由引理 1.3，二者基频分别为 $\tfrac{2}{\pi}\sin(2\pi ft)$ 与 $\tfrac{2}{\pi}\sin(2\pi ft+\tfrac\pi2)=\tfrac{2}{\pi}\cos(2\pi ft)$。参考取 $r_I=\sin(2\pi ft)$、$r_Q=\cos(2\pi ft)$。

### 定理 5.1（I/Q 零串扰分离）

$$
\mathrm{LPF}[s\,r_I]=\frac{w_I}{\pi}R_I+\mathcal{O}\!\big(\text{阻带}\big),\qquad
\mathrm{LPF}[s\,r_Q]=\frac{w_Q}{\pi}R_Q+\mathcal{O}\!\big(\text{阻带}\big),
$$

即同频两路互不串扰，重建 $\hat R_I=\pi\,\mathrm{LPF}[s r_I]/w_I$，$\hat R_Q=\pi\,\mathrm{LPF}[s r_Q]/w_Q$。

**证明.** 仅看该频率组（其余频率由引理 2.1 归零）。
- $r_I$ 上：I 基频 $\langle\sin\cdot\sin\rangle=\tfrac12$ → $\tfrac{w_I}{\pi}R_I$；Q 基频 $\langle\cos\cdot\sin\rangle=0$（推论 2.2）→ 无 Q 泄漏。
- $r_Q$ 上：对称地得 $\tfrac{w_Q}{\pi}R_Q$，I 泄漏为 $\langle\sin\cdot\cos\rangle=0$。
- **奇次谐波不在 $f$ 处**：I、Q 的谐波位于 $3f,5f,\dots$，与参考频率 $f$ 不等，引理 2.1 归零；偶次谐波已被 $d=\tfrac12$ 消去（引理 1.3）。

故同频两路在基频上严格正交，无串扰。$\blacksquare$

### 5.2 编码矩阵与条件数

将一个 IQ 对的"基频投影"写成线性映射：测量向量 $\mathbf{y}=(\,\mathrm{LPF}[s r_I],\ \mathrm{LPF}[s r_Q]\,)^\top$ 与反射率 $\mathbf{R}=(R_I,R_Q)^\top$ 满足

$$
\mathbf{y}=\mathbf{A}\,\mathbf{R},\qquad
\mathbf{A}=\frac{1}{\pi}\begin{pmatrix}w_I & 0\\ 0 & w_Q\end{pmatrix}.
$$

由定理 5.1 的正交性，非对角元为零。其（白化后，即以 $w$ 归一）条件数

$$
\kappa(\mathbf{A})=\frac{\sigma_{\max}}{\sigma_{\min}}=1.
$$

**意义.** $\kappa=1$ 表示编码矩阵完美正交，求逆不放大噪声——这是 $d=\tfrac12$ 相位编码相对占空比编码（$d=0.25/0.75$，存在偶次谐波、$\kappa>1$）的根本优势。

### 5.3 串扰来源小结

唯一的同频串扰来自：(i) 慢变假设残差 $O(B_R/f)$（引理 2.3）；(ii) 积分时延导致的 I↔Q 旋转——这正是 §6 必须补偿的原因。

---

## §6 传感器积分时延与相位补偿

### 引理 6.1（boxcar 积分的群延迟与幅值衰减）

对单频载波 $x(t)=\sin(2\pi f t+\psi)$，前向 boxcar 算子 (3.2) 给出

$$
(\mathcal{S}x)(t)=\mathrm{sinc}(\pi f T_\mathrm{int})\,\sin\!\big(2\pi f(t+\tau)+\psi\big),\qquad
\boxed{\;\tau=\frac{T_\mathrm{int}}{2}\;},
$$

其中 $\mathrm{sinc}(u)=\sin u/u$。即积分使载波**相位超前** $\Delta\phi=2\pi f\tau=\pi f T_\mathrm{int}$，**幅值衰减** $\mathrm{sinc}(\pi f T_\mathrm{int})$。

**证明.**

$$
(\mathcal{S}x)(t)=\frac{1}{T_\mathrm{int}}\int_t^{t+T_\mathrm{int}}\sin(2\pi f\tau'+\psi)\,d\tau'
=\frac{-1}{2\pi f T_\mathrm{int}}\Big[\cos(2\pi f\tau'+\psi)\Big]_t^{t+T_\mathrm{int}}.
$$

用 $\cos\beta-\cos\alpha=-2\sin\tfrac{\alpha+\beta}{2}\sin\tfrac{\beta-\alpha}{2}$，令 $\alpha=2\pi ft+\psi$，$\beta=\alpha+2\pi fT_\mathrm{int}$：

$$
=\frac{1}{2\pi f T_\mathrm{int}}\cdot 2\sin\!\Big(2\pi f t+\psi+\pi f T_\mathrm{int}\Big)\sin\!\big(\pi f T_\mathrm{int}\big)
=\mathrm{sinc}(\pi f T_\mathrm{int})\,\sin\!\big(2\pi f(t+\tfrac{T_\mathrm{int}}{2})+\psi\big).
$$

即得。**离散修正**：对 $n_\mathrm{per}$ 点等间隔求和，群延迟为质心 $\tau=\dfrac{(n_\mathrm{per}-1)}{2F_\mathrm{sim}}$（代码所用），当 $n_\mathrm{per}$ 大时 $\to T_\mathrm{int}/2$。$\blacksquare$

### 定理 6.2（未补偿时延的幅值损失与 I/Q 旋转）

若用未超前的参考 $r_I=\sin(2\pi ft)$ 解调被延迟的 I 通道载波 $\propto\sin(2\pi f(t+\tau))$，则

$$
\mathrm{LPF}[\,\cdot\,]\ \propto\ \tfrac12\big[\underbrace{\cos\Delta\phi}_{\text{保留}}\,R_I\ +\ \underbrace{\sin\Delta\phi}_{\text{泄漏入 }Q}\,(\text{正交分量})\big],\quad \Delta\phi=2\pi f\tau.
$$

即：幅值按 $\cos\Delta\phi$ 衰减，并有 $\sin\Delta\phi$ 比例的能量旋转进入正交（Q）通道。

**证明.** 由引理 2.1，$\langle\sin(2\pi f(t+\tau))\sin(2\pi ft)\rangle=\tfrac12\cos(2\pi f\tau)=\tfrac12\cos\Delta\phi$；而与 $\cos(2\pi ft)$ 的内积为 $\tfrac12\sin\Delta\phi$。前者为同相保留，后者为正交泄漏。$\blacksquare$

> **数值（本方案，$\tau\approx2.4$–$2.5$ ms）**：$f=53$ Hz 时 $\Delta\phi=2\pi\cdot53\cdot0.0024\approx0.80\ \mathrm{rad}=45.6^\circ$，未补偿则幅值仅剩 $\cos45.6^\circ\approx70\%$，且 $\sin45.6^\circ\approx71\%$ 的能量串入正交通道——高频通道几乎失效。

### 定理 6.3（时延补偿的最优性）

将参考**超前** $\tau$，即取 $r_m(t)=\sin\!\big(2\pi f_m(t+\tau)+\varphi_m\big)$，则解调中等效相位差 $\Delta\phi\to0$，定理 4.1 与定理 5.1 的结论精确成立（同相幅值 $\cos0=1$，正交泄漏 $\sin0=0$）。在所有"参考时移"中，超前量 $\tau=T_\mathrm{int}/2$ 是唯一使同相幅值取极大且正交泄漏为零者。

**证明.** 由引理 6.1，被测载波相位为 $2\pi f(t+\tau)+\psi$。令参考相位与之一致即 $r_m=\sin(2\pi f(t+\tau)+\varphi_m)$，则定理 6.2 中的等效 $\Delta\phi=0$，$\cos\Delta\phi$ 取唯一极大 $1$、$\sin\Delta\phi=0$。$\blacksquare$

代码实现（`spectral_reconstruction.py` 的 `_reference` 与 `group_delay`）：

```python
tau = integration_time / 2          # 或显式 integration_delay
ref = sin(2*pi*f*(t + tau) + phase) # 超前补偿后的匹配参考
```

> **幅值项 $\mathrm{sinc}(\pi f T_\mathrm{int})$** 不随相位补偿消失，但它是各通道已知的确定衰减，被 $w_k$ 校准或参考标定（§9）一并吸收，不影响重建。

---

## §7 多频串扰与频率规划约束

定理 4.1 第 4 步要求跨频残差落入 $\mathrm{LPF}$ 阻带。下面给出可操作的设计约束。

### 命题 7.1（基频拍频隔离）

用 $r_i$（频率 $f_i$）解调他通道 $j$（频率 $f_j$）的基频时，乘积出现于 $|f_i\pm f_j|$。差频分量 $|f_i-f_j|$ 必须落在阻带：

$$
\boxed{\,|f_i-f_j|>f_c\,}\quad(\text{推荐留余量 }\ge 2f_c).
$$

否则 $j$ 的能量以拍频 $|f_i-f_j|<f_c$ 漏入通道 $i$。

### 命题 7.2（谐波混叠隔离）

$d=\tfrac12$ 仅余奇次谐波 $kf_i\ (k=3,5,\dots)$。采样率 $F_s$ 下它们混叠到

$$
\mathrm{alias}(kf_i)=\big|\,kf_i-\mathrm{round}(kf_i/F_s)\,F_s\,\big|.
$$

要求对所有奇 $k$ 与所有 $j$：

$$
\boxed{\,\big|\mathrm{alias}(kf_i)-f_j\big|>f_c\,.}
$$

### 命题 7.3（奈奎斯特约束）

$$
\boxed{\,f_{\max}=\max_m f_m<F_s/2\,.}
$$

### 7.4 本方案取值的相容性

频率组 $\{13,23,33,43,53\}$ Hz（间距 $G=10$ Hz）满足：

- 基频隔离：相邻间距 $10>2f_c=7$ ✔
- 奈奎斯特：$f_{\max}=53<100$ ✔，且 §6 相位误差 $45.6^\circ$ 可补偿 ✔
- 谐波：如 $13\times3=39$ Hz，距 $33$ Hz 为 $6$ Hz（$>f_c=3.5$，余量偏小，构成 628/717 nm 通道的主要系统误差源——见 §10.3）。

> **设计权衡.** 增大 $G$ 抬高 $f_{\max}$（相位误差与 $\mathrm{sinc}$ 衰减加剧），减小 $G$ 触犯基频/谐波隔离。实验扫描（系统文档 §4.3）证明 $G=10$ Hz 为综合 RMSE 最小点。

---

## §8 时间分辨率

### 命题 8.1

锁相系统可无失真追踪的反射率最高频率等于 $\mathrm{LPF}$ 截止 $f_c$（4 阶 Butterworth 在 $f_c$ 处约 $-3$ dB）；其零相位由双向滤波 `sosfiltfilt` 保证。同时必须满足

$$
B_R\le f_c<\tfrac12\min_{i\neq j}|f_i-f_j|,
$$

即追踪带宽不得超过频率间隔之半，否则相邻通道的拍频边带进入通带造成串扰。本方案 $f_c=3.5$ Hz，目标最高变化频率 $2.55$ Hz，$f_c/2.55\approx1.37$ 倍余量；而 $\tfrac12\min|f_i-f_j|=5$ Hz $>f_c$，约束相容。

---

## §9 光谱加权、波长重采样与参考标定

### 9.1 有效权重

第 $k$ 通道的 LED–传感器有效权重为发射光谱与响应的重叠积分

$$
w_k=\int_{\Lambda} I_k(\lambda)\,\eta(\lambda)\,d\lambda.\tag{9.1}
$$

代码默认 $I_k$ 为高斯线型 $I_k(\lambda)=\exp\!\big[-4\ln2\,(\lambda-\lambda_k)^2/\Delta\lambda_k^2\big]$（FWHM $\Delta\lambda_k$）。

### 9.2 波长重采样（处理采样间隔不一致）

实际中 $I_k$（如数据手册，间隔 5 nm）与 $\eta$（如标定曲线，间隔 1 或 10 nm）常在**不同波长栅格**上给出。设二者支撑并集为 $[\lambda_{\min},\lambda_{\max}]$，取统一步长 $\delta$（`resample_step`，默认 1 nm）的细栅格 $\Lambda_\delta$，以线性插值（栅格外补零）将两者重采样到 $\Lambda_\delta$，再用梯形法计算 (9.1)：

$$
w_k\approx\sum_{\lambda\in\Lambda_\delta}\tilde I_k(\lambda)\,\tilde\eta(\lambda)\,\delta.
$$

实现见 `overlap_integral()`（`code/spectral_reconstruction.py`）。该步骤保证不同源数据可一致积分。

### 9.3 标定模式（按实用性排序）

$\eta(\lambda)$ 常**不可得**，故权重 $w_k$ 未必已知。系统支持四种标定：

| 模式 | 条件 | 输出 |
|------|------|------|
| `white_reference` | 拍一帧已知反射率 $R_\mathrm{ref}$ 的平整参考板 | **绝对**反射率 |
| `weights` | 已标定的 $w_k$ | 绝对反射率 |
| `spectral` | 已知 $\eta(\lambda)$（可选） | **相对**反射率 |
| `none` | 三者皆无 | 至未知逐通道尺度的相对量 |

### 定理 9.4（参考板标定消去未知量）

设参考板反射率为已知常数 $R_\mathrm{ref}$（白板 $R_\mathrm{ref}=1$，或灰板 $R_\mathrm{ref}=0.5$）。其锁相输出为

$$
A_{\mathrm{ref},k}=\mathrm{LPF}[\,s_\mathrm{ref}\,r_k\,]=\frac{w_k}{\pi}R_\mathrm{ref}\cdot \underbrace{\mathrm{sinc}(\pi f_k T_\mathrm{int})\cos\Delta\phi_k}_{\text{已知确定衰减}}.
$$

对测量信号同理 $y_k=\tfrac{w_k}{\pi}R_k\cdot(\text{相同衰减})$。则

$$
\boxed{\;\hat R_k(t)=R_\mathrm{ref}\,\frac{y_k(t)}{A_{\mathrm{ref},k}}\;}
$$

**精确消去** $w_k$、常数 $\pi$、$\mathrm{sinc}$ 幅值衰减乃至残余相位因子 $\cos\Delta\phi_k$——无需知道 LED 功率与传感器响应。

**证明.** 两式相除，所有与 $t$ 无关的通道因子（$w_k,\pi,\mathrm{sinc},\cos\Delta\phi_k$）成比例抵消，仅余 $R_k(t)/R_\mathrm{ref}$，乘 $R_\mathrm{ref}$ 即得。$\blacksquare$

> **工程提示.** $A_{\mathrm{ref},k}$ 应取参考板锁相输出在内区的稳健标量（如中位数），而非逐点时间序列相除——后者会因 IIR 边缘振铃经过零点而发散（实现中已用内区中位数）。参考板须使各 LED 同亮且**不饱和 ADC**：若所有通道 $R=1$ 同时点亮可能超满量程，建议用灰板（如 $R_\mathrm{ref}=0.5$）。

---

## §10 误差分析

### 10.1 锁相对噪声的抑制

设 $\varepsilon[n]$ 为方差 $\sigma_\varepsilon^2$ 的白噪声。混频 $\varepsilon\cdot r_k$（$|r_k|\le1$）不改变其白谱密度量级；经等效噪声带宽 $B_n$ 的 $\mathrm{LPF}$（对 4 阶 Butterworth，$B_n\approx1.026\,f_c$）后，输出噪声方差

$$
\sigma_y^2\approx\sigma_\varepsilon^2\,\frac{B_n}{F_s/2}\cdot\frac12 .
$$

反射率噪声经定理 4.1 放大 $\pi/w_k$：

$$
\boxed{\;\sigma_{\hat R_k}\approx\frac{\pi}{w_k}\,\sigma_\varepsilon\sqrt{\frac{B_n}{F_s}}\;}
$$

**意义.** 噪声随 $\sqrt{f_c/F_s}$ 下降：窄带 $\mathrm{LPF}$（小 $f_c$）与高采样率 $F_s$ 均抑噪，但 $f_c$ 受时间分辨率（§8）下限约束，二者需折中。权重 $w_k$ 越小（如 450、850 nm 边缘通道）噪声放大越大。

### 10.2 量化误差

$B$ 位均匀量化等效附加均匀噪声，方差 $\sigma_q^2=q^2/12=2^{-2B}/12$。本方案 $B=8$：$\sigma_q=2^{-8}/\sqrt{12}\approx1.1\times10^{-3}$，小于电子噪声 $\sigma_\varepsilon\approx2.1\times10^{-3}$（SNR $=42$ dB），故**量化非瓶颈**，8 位足够。

### 10.3 系统误差（偏差）

非随机误差主要来自：

1. **谐波混叠残余**（命题 7.2）：如 $53\times3=159$ Hz 混叠近 $43$ Hz（差 $2$ Hz $<f_c$），在 628/717 nm 通道引入 $\mathrm{RMSE}\approx0.05$–$0.06$ 的系统偏差。
2. **慢变假设残差**（引理 2.3）：$O(B_R/f_m)$，对高 $f_2$（接近 $f_c$）的快变通道（717–850 nm）更明显。
3. **IIR 边缘瞬态**：双向 Butterworth 在信号首尾产生约 $1$ s 过渡段，需 `trim`。

### 10.4 与实测一致性

按上述模型，10 通道均值 $\mathrm{RMSE}=0.044$、皮尔逊 $r=0.964$。本仓库 `tests/test_reconstruction.py` 用前向仿真驱动配置化 API，复现**逐通道**数值（450 nm $0.0246$、850 nm $0.0717$、均值 $0.0439$），与理论及原始实现一致。

---

## §11 定理汇总与设计准则

**重建主公式**（定理 4.1 / 5.1 / 6.3 / 9.4）：

$$
\hat R_k(t)=\frac{\pi}{w_k}\,\mathrm{LPF}\Big[s(t)\,\sin\!\big(2\pi f_k(t+\tau)+\varphi_k\big)\Big],
\qquad\tau=\frac{T_\mathrm{int}}{2},
$$

或免标定的参考板形式 $\hat R_k=R_\mathrm{ref}\,y_k/A_{\mathrm{ref},k}$。

**设计准则**：

| 准则 | 来源 | 本方案取值 |
|------|------|-----------|
| $d=\tfrac12$（消偶次谐波、$\kappa=1$、基频最大） | 引理 1.3、§5.2 | $0.5$ |
| I/Q 相差 $90^\circ$ | 推论 2.2、定理 5.1 | $0^\circ,90^\circ$ |
| 参考超前 $\tau=T_\mathrm{int}/2$ | 定理 6.3 | $2.4$–$2.5$ ms |
| 基频隔离 $|f_i-f_j|>f_c$（余量 $2f_c$） | 命题 7.1 | 间距 $10>7$ Hz |
| 谐波隔离 $|\mathrm{alias}(kf_i)-f_j|>f_c$ | 命题 7.2 | 经扫描优化 |
| 奈奎斯特 $f_{\max}<F_s/2$ | 命题 7.3 | $53<100$ Hz |
| 追踪带宽 $B_R\le f_c<\tfrac12\min|f_i-f_j|$ | 命题 8.1 | $2.55\le3.5<5$ Hz |
| 位深 $B$ 使 $\sigma_q\ll\sigma_\varepsilon$ | §10.2 | $8$ bit |

---

## 参考文献

1. R. N. Bracewell, *The Fourier Transform and Its Applications*, 3rd ed., McGraw-Hill, 2000.（傅里叶级数、卷积定理）
2. A. V. Oppenheim, R. W. Schafer, *Discrete-Time Signal Processing*, 3rd ed., Prentice Hall, 2009.（采样、混叠、boxcar 群延迟）
3. P. A. Temple, "An introduction to phase-sensitive amplifiers (lock-in detection)," *Am. J. Phys.* 43(9), 1975.（锁相检测原理）
4. J. H. Scofield, "A frequency-domain description of a lock-in amplifier," *Am. J. Phys.* 62(2), 1994.（锁相噪声带宽）
5. W. Press et al., *Numerical Recipes*, 3rd ed., 2007.（量化噪声、数值积分）

---

*文档结束 — 配套代码与测试见本仓库 `code/` 与 `tests/`。*
