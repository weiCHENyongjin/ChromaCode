"""Forward simulator for the IQ phase-encoded multispectral sensing system.

多波长 LED IQ 相位编码光强采集系统 —— 前向仿真器。

This module generates a synthetic single-detector capture for the validated v1
configuration, then reconstructs it through :class:`SpectralReconstructor`
(the shared inverse in ``spectral_reconstruction.py``) so the forward and
inverse stay in lock-step with a single demodulation implementation.

Validated configuration (v1)
----------------------------
* Sensor: ``fs = 200 Hz``, 8-bit, SNR ≈ 42 dB.
* LEDs: ``[13, 23, 33, 43, 53] Hz`` × duty 0.5 × phases ``{0°, 90°}`` = 10 channels.
* Wavelengths: ``[450, 494, 539, 583, 628, 672, 717, 761, 806, 850] nm``.
* Demodulation: 4th-order Butterworth LPF (3.5 Hz) + integration-delay
  compensation.
* Performance: mean RMSE ≈ 0.044, mean Pearson r ≈ 0.964.

Run ``python code/iq_sensing_system.py`` to reproduce the result and figure.
See ``../CODE_STYLE.md`` for the conventions used throughout.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, replace
from typing import Dict, List

import numpy as np
import scipy.signal as signal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spectral_reconstruction import (  # noqa: E402  (local import after path setup)
    DEFAULT_FWHM_NM,
    DEFAULT_LPF_CUTOFF_HZ,
    LEDChannel,
    ReconstructionConfig,
    ReconstructionResult,
    SpectralReconstructor,
    gaussian_emission,
    silicon_sensor_response,
)


# ─────────────────────────────────────────────────────────────────────────
# Constants: the validated 10-channel layout
# ─────────────────────────────────────────────────────────────────────────
#: Channel center wavelengths (nm), in channel order.
WAVELENGTHS_NM: List[float] = [450, 494, 539, 583, 628, 672, 717, 761, 806, 850]

#: LED modulation frequencies (Hz); each carries an I and a Q channel.
MODULATION_FREQUENCIES_HZ: List[float] = [13, 23, 33, 43, 53]

#: Initial phases (radians) for the I (sin) and Q (cos) channels of each group.
PHASES_RAD = (0.0, np.pi / 2)

#: Spectral grid (nm) used for weight integrals.
SPECTRAL_GRID_NM: np.ndarray = np.arange(400, 1001, 1.0)


# ─────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────
@dataclass
class SimulationParameters:
    """Numeric parameters of the forward simulation.

    Attributes
    ----------
    sample_rate_hz : float
        Detector sample rate ``fs``.
    oversample_rate_hz : float
        High-rate grid used to synthesize the continuous-time signal
        (simulation only; real hardware needs no oversampling).
    duration_s : float
        Total simulated duration.
    snr_db : float
        Electronic-noise signal-to-noise ratio.
    adc_bits : int
        ADC bit depth.
    lpf_cutoff_hz : float
        Lock-in low-pass cutoff used for the inverse.
    """
    sample_rate_hz: float = 200.0
    oversample_rate_hz: float = 5000.0
    duration_s: float = 12.0
    snr_db: float = 42.0
    adc_bits: int = 8
    lpf_cutoff_hz: float = DEFAULT_LPF_CUTOFF_HZ

    @property
    def samples_per_period(self) -> int:
        """Oversampled points averaged per detector sample (boxcar width)."""
        return int(self.oversample_rate_hz // self.sample_rate_hz)

    @property
    def n_samples(self) -> int:
        """Number of detector samples."""
        return int(self.duration_s * self.sample_rate_hz)

    @property
    def adc_max(self) -> int:
        """Largest ADC code (full scale)."""
        return 2 ** self.adc_bits - 1

    @property
    def boxcar_delay_s(self) -> float:
        """Group delay (s) of the discrete boxcar integration; see theory §6."""
        return (self.samples_per_period - 1) / (2 * self.oversample_rate_hz)


# ─────────────────────────────────────────────────────────────────────────
# System construction
# ─────────────────────────────────────────────────────────────────────────
def default_channels() -> List[LEDChannel]:
    """Build the validated 10-channel LED layout (5 frequencies × IQ phases)."""
    channels: List[LEDChannel] = []
    for index, wavelength_nm in enumerate(WAVELENGTHS_NM):
        frequency_hz = MODULATION_FREQUENCIES_HZ[index // 2]
        phase_rad = PHASES_RAD[index % 2]
        channels.append(LEDChannel(
            wavelength_nm=float(wavelength_nm),
            frequency_hz=float(frequency_hz),
            phase_rad=phase_rad,
            fwhm_nm=DEFAULT_FWHM_NM,
            label=f"ch{index + 1:02d}",
        ))
    return channels


def channel_weights(channels: List[LEDChannel],
                    scale: float = 0.2) -> Dict[float, float]:
    """LED-to-sensor effective weights, normalized so the largest equals ``scale``.

    Each weight is the spectral overlap of the LED line with the silicon sensor
    response. ``scale`` sets the largest channel's share of full scale.

    Parameters
    ----------
    channels : list of LEDChannel
    scale : float
        Peak weight after normalization (default 0.2 ≈ 20% of full scale).

    Returns
    -------
    dict
        Maps ``wavelength_nm`` to its effective weight.
    """
    response = silicon_sensor_response(SPECTRAL_GRID_NM)
    raw = np.array([
        np.trapezoid(gaussian_emission(SPECTRAL_GRID_NM, ch.wavelength_nm,
                                       ch.fwhm_nm) * response, SPECTRAL_GRID_NM)
        for ch in channels
    ])
    scaled = raw / raw.max() * scale
    return {ch.wavelength_nm: float(scaled[i]) for i, ch in enumerate(channels)}


def synthetic_target_reflectance(time_s: np.ndarray,
                                 wavelength_nm: float) -> np.ndarray:
    """Synthetic time-varying reflectance of the target at one wavelength.

    The baseline varies with wavelength (Gaussian peak near 600 nm); two sine
    components model dynamics. The change frequencies increase with wavelength:
    ``f1 = 0.5 + (λ − 450)/400`` Hz and ``f2 = 1.7 f1``.

    Parameters
    ----------
    time_s : np.ndarray
        Time axis (s).
    wavelength_nm : float
        Wavelength (nm).

    Returns
    -------
    np.ndarray
        Reflectance in [0.05, 0.95].
    """
    f1_hz = 0.5 + (wavelength_nm - 450) / 400
    f2_hz = 1.7 * f1_hz
    reflectance = (0.4
                   + 0.3 * np.exp(-((wavelength_nm - 600) / 150) ** 2)
                   + 0.15 * np.sin(2 * np.pi * f1_hz * time_s)
                   + 0.08 * np.sin(2 * np.pi * f2_hz * time_s + np.pi / 3))
    return np.clip(reflectance, 0.05, 0.95)


# ─────────────────────────────────────────────────────────────────────────
# Forward simulator
# ─────────────────────────────────────────────────────────────────────────
class ForwardSimulator:
    """Synthesize a single-detector capture for a multi-LED system."""

    def __init__(self, channels: List[LEDChannel],
                 params: SimulationParameters,
                 weight_scale: float = 0.2) -> None:
        self.channels = channels
        self.params = params
        self.weights = channel_weights(channels, weight_scale)

    @classmethod
    def default(cls) -> "ForwardSimulator":
        """The validated v1 system (10 channels, default parameters)."""
        return cls(default_channels(), SimulationParameters())

    # ----- time axes -------------------------------------------------------
    def oversample_time_s(self) -> np.ndarray:
        """High-rate time axis used to synthesize the continuous signal."""
        n = int(self.params.duration_s * self.params.oversample_rate_hz)
        return np.linspace(0, self.params.duration_s, n, endpoint=False)

    def sensor_time_s(self) -> np.ndarray:
        """Detector sample time axis."""
        return np.linspace(0, self.params.duration_s,
                           self.params.n_samples, endpoint=False)

    # ----- forward model ---------------------------------------------------
    def synthesize(self,
                   reflectance_by_wavelength: Dict[float, np.ndarray],
                   seed: int = 42) -> np.ndarray:
        """Synthesize the quantized detector signal.

        Steps: (1) sum each LED's PWM contribution on the oversampled grid,
        (2) boxcar-integrate to the detector rate, (3) add Gaussian electronic
        noise at ``snr_db``, (4) quantize to ``adc_bits``.

        Parameters
        ----------
        reflectance_by_wavelength : dict
            Maps ``wavelength_nm`` to the true reflectance on the oversampled
            time axis (see :meth:`oversample_time_s`).
        seed : int
            RNG seed for the additive noise (reproducibility).

        Returns
        -------
        np.ndarray
            The quantized detector signal on :meth:`sensor_time_s`.
        """
        params = self.params
        # Legacy global-seed RNG is kept intentionally: it reproduces the exact
        # per-channel RMSE values published in the docs/README. 保持发布数值一致。
        np.random.seed(seed)
        time_hi = self.oversample_time_s()
        n_samples = params.n_samples
        samples_per_period = params.samples_per_period

        # (1) Continuous-time superposition of every LED's PWM contribution.
        continuous = np.zeros(time_hi.size)
        for channel in self.channels:
            pwm = (signal.square(2 * np.pi * channel.frequency_hz * time_hi
                                 + channel.phase_rad, duty=0.5) + 1) / 2
            continuous += (self.weights[channel.wavelength_nm] * pwm
                           * reflectance_by_wavelength[channel.wavelength_nm])

        # (2) Boxcar integration: average samples_per_period points per sample.
        usable = n_samples * samples_per_period
        sampled = continuous[:usable].reshape(n_samples, samples_per_period).mean(axis=1)

        # (3) Additive Gaussian electronic noise at the requested SNR.
        rms = np.sqrt(np.mean(sampled ** 2))
        noise_std = rms / 10 ** (params.snr_db / 20)
        noisy = sampled + np.random.normal(0, noise_std, n_samples)

        # (4) Uniform ADC quantization.
        adc_max = params.adc_max
        return np.clip(np.round(noisy * adc_max), 0, adc_max) / adc_max

    # ----- matching inverse ------------------------------------------------
    def reconstruction_config(self) -> ReconstructionConfig:
        """A config that inverts this simulator exactly (known weights + delay)."""
        calibrated = [replace(ch, weight=self.weights[ch.wavelength_nm])
                      for ch in self.channels]
        return ReconstructionConfig(
            sample_rate_hz=self.params.sample_rate_hz,
            channels=calibrated,
            lpf_cutoff_hz=self.params.lpf_cutoff_hz,
            integration_delay_s=self.params.boxcar_delay_s,
        )

    def reconstructor(self) -> SpectralReconstructor:
        """A :class:`SpectralReconstructor` matched to this simulator."""
        return SpectralReconstructor(self.reconstruction_config())


# ─────────────────────────────────────────────────────────────────────────
# Evaluation
# ─────────────────────────────────────────────────────────────────────────
def evaluate_reconstruction(result: ReconstructionResult,
                            oversample_time_s: np.ndarray,
                            true_reflectance: Dict[float, np.ndarray],
                            ) -> Dict[float, Dict[str, float]]:
    """RMSE and Pearson r per channel against the ground truth.

    Parameters
    ----------
    result : ReconstructionResult
        The reconstruction to score.
    oversample_time_s : np.ndarray
        Time axis on which ``true_reflectance`` is defined.
    true_reflectance : dict
        Ground-truth reflectance by wavelength on ``oversample_time_s``.

    Returns
    -------
    dict
        Maps ``wavelength_nm`` to ``{"rmse": ..., "corr": ...}``.
    """
    metrics: Dict[float, Dict[str, float]] = {}
    for wavelength_nm, series in result.reflectance.items():
        truth = np.interp(result.time_s, oversample_time_s,
                          true_reflectance[wavelength_nm])
        rmse = float(np.sqrt(np.nanmean((series - truth) ** 2)))
        corr = float(np.corrcoef(series, truth)[0, 1])
        metrics[wavelength_nm] = {"rmse": rmse, "corr": corr}
    return metrics


# ─────────────────────────────────────────────────────────────────────────
# Demo / CLI
# ─────────────────────────────────────────────────────────────────────────
def _plot(result: ReconstructionResult, simulator: ForwardSimulator,
          time_hi: np.ndarray, true_reflectance: Dict[float, np.ndarray],
          metrics: Dict[float, Dict[str, float]], out_path: str) -> None:
    """Save a 5×2 grid of ground-truth vs reconstructed reflectance.

    English labels keep this portable on systems without a CJK font; for the
    polished bilingual README figures use ``code/make_figures.py``.
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(5, 2, figsize=(14, 18), sharex=True)
    for index, channel in enumerate(simulator.channels):
        ax = axes[index // 2, index % 2]
        wavelength_nm = channel.wavelength_nm
        truth = np.interp(result.time_s, time_hi, true_reflectance[wavelength_nm])
        phase = "I" if channel.phase_rad < 0.1 else "Q"
        ax.plot(result.time_s, truth, lw=2, alpha=0.6, label="ground truth")
        ax.plot(result.time_s, result.reflectance[wavelength_nm], "k--", lw=1.5,
                label=f"reconstructed RMSE={metrics[wavelength_nm]['rmse']:.4f}")
        ax.set_ylim(0, 1.1)
        ax.set_title(f"{wavelength_nm:.0f} nm | {channel.frequency_hz:.0f} Hz-{phase}",
                     fontsize=10)
        ax.legend(fontsize=8)
    axes[-1, 0].set_xlabel("time (s)")
    axes[-1, 1].set_xlabel("time (s)")
    fig.suptitle("IQ phase-encoded system - 10-channel reconstruction",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Reproduce the validated reconstruction and save a result figure."""
    simulator = ForwardSimulator.default()
    time_hi = simulator.oversample_time_s()
    true_reflectance = {wl: synthetic_target_reflectance(time_hi, wl)
                        for wl in WAVELENGTHS_NM}

    print("synthesizing detector signal ...")
    signal_q = simulator.synthesize(true_reflectance)
    sensor_time = simulator.sensor_time_s()

    print("IQ lock-in reconstruction ...")
    result = simulator.reconstructor().reconstruct(
        signal_q, time_s=sensor_time, trim_s=1.0)
    metrics = evaluate_reconstruction(result, time_hi, true_reflectance)

    print(f"\n{'channel':<10}{'freq':>8}{'phase':>7}{'RMSE':>11}{'Pearson r':>11}")
    print("-" * 47)
    for index, channel in enumerate(simulator.channels):
        phase = "I" if channel.phase_rad < 0.1 else "Q"
        m = metrics[channel.wavelength_nm]
        print(f"{channel.wavelength_nm:>5.0f} nm   {channel.frequency_hz:>5.0f}Hz"
              f"{phase:>7}{m['rmse']:>11.5f}{m['corr']:>11.5f}")
    mean_rmse = np.mean([m["rmse"] for m in metrics.values()])
    mean_corr = np.mean([m["corr"] for m in metrics.values()])
    print("-" * 47)
    print(f"{'mean':<10}{'':>8}{'':>7}{mean_rmse:>11.5f}{mean_corr:>11.5f}")

    try:
        out_path = "reconstruction_result.png"
        _plot(result, simulator, time_hi, true_reflectance, metrics, out_path)
        print(f"\nresult figure saved: {out_path}")
    except Exception as exc:  # pragma: no cover - plotting is optional
        print(f"(plot skipped: {exc})")


if __name__ == "__main__":
    main()
