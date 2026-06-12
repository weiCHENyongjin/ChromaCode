"""Config-driven multispectral reconstruction via IQ phase-encoded lock-in.

配置驱动的多光谱重建（IQ 相位编码锁相解调）—— 本项目的核心库。

This module is the **single source of truth** for the inverse problem and the
spectral model. The forward simulator (``iq_sensing_system.py``) reuses the
classes defined here rather than re-implementing the demodulation.

Pipeline
--------
1. Describe the **sensor** and the **LEDs** in a JSON/YAML config
   (:class:`ReconstructionConfig`).
2. Feed the captured 1-D detector time series to
   :meth:`SpectralReconstructor.reconstruct`.
3. Receive the reconstructed time-varying reflectance for every wavelength
   (:class:`ReconstructionResult`).

Calibration
-----------
The sensor spectral response is usually **unknown**, so four calibration modes
are supported (in order of practicality):

``white_reference``
    Supply one capture of a flat reference target of known reflectance. Yields
    absolute reflectance with no knowledge of LED power or sensor response.
``weights``
    Per-channel weights known from a prior calibration.
``spectral``
    Compute weights from the (optional) sensor spectral response and the LED
    spectra. Yields *relative* reflectance (absolute scale unknown).
``none``
    No calibration data: reflectance up to an unknown per-channel scale.

Mismatched wavelength sampling between LED spectra and the response curve is
handled by resampling onto a common grid (see :func:`spectral_overlap`).

Dependencies
------------
``numpy``, ``scipy``; YAML configs additionally require ``pyyaml``.

See ``../CODE_STYLE.md`` for the conventions used throughout.
"""

from __future__ import annotations

import json
import os
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.signal import butter, sosfiltfilt

try:
    import yaml  # optional; only needed for .yaml/.yml configs
    _HAS_YAML = True
except Exception:  # pragma: no cover - environment dependent
    _HAS_YAML = False


# ─────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────
#: Default common-grid step (nm) for spectral overlap integrals.
DEFAULT_RESAMPLE_STEP_NM: float = 1.0

#: Default lock-in low-pass cutoff (Hz) and Butterworth order.
DEFAULT_LPF_CUTOFF_HZ: float = 3.5
DEFAULT_LPF_ORDER: int = 4

#: Default LED line-width (nm), full width at half maximum.
DEFAULT_FWHM_NM: float = 20.0


# ─────────────────────────────────────────────────────────────────────────
# Spectral model (shared with the forward simulator)
# ─────────────────────────────────────────────────────────────────────────
def gaussian_emission(grid_nm: np.ndarray,
                      center_nm: float,
                      fwhm_nm: float) -> np.ndarray:
    """Gaussian emission line with unit peak.

    Parameters
    ----------
    grid_nm : np.ndarray
        Wavelength axis (nm) on which to evaluate the line.
    center_nm : float
        Center wavelength (nm).
    fwhm_nm : float
        Full width at half maximum (nm).

    Returns
    -------
    np.ndarray
        Line shape evaluated on ``grid_nm``, peaking at 1.0.
    """
    return np.exp(-4.0 * np.log(2.0) * ((grid_nm - center_nm) / fwhm_nm) ** 2)


def silicon_sensor_response(grid_nm: np.ndarray) -> np.ndarray:
    """Approximate, normalized spectral response of a silicon detector.

    A smooth trapezoid-like curve over 400-1000 nm, peak-normalized to 1.0.
    Used as a stand-in when a measured response is unavailable.

    Parameters
    ----------
    grid_nm : np.ndarray
        Wavelength axis (nm).

    Returns
    -------
    np.ndarray
        Normalized response (0..1) on ``grid_nm``; zero outside 400-1000 nm.
    """
    response = np.zeros_like(grid_nm, dtype=float)
    in_band = (grid_nm >= 400) & (grid_nm <= 1000)
    wl = grid_nm[in_band]
    # Rising edge to ~500 nm, falling edge past ~700 nm (silicon-like).
    response[in_band] = (np.clip((wl - 380) / 120, 0, 1)
                         * np.clip((1020 - wl) / 300, 0, 1))
    peak = response.max()
    return response / peak if peak > 0 else response


def spectral_overlap(channel: "LEDChannel",
                     response_wavelengths_nm: np.ndarray,
                     response_values: np.ndarray,
                     resample_step_nm: float = DEFAULT_RESAMPLE_STEP_NM) -> float:
    """Spectral overlap integral of one LED with the sensor response.

    Both curves are resampled onto a common fine grid before integration, so
    their native wavelength sampling intervals may differ (e.g. a 5 nm
    datasheet spectrum vs a 10 nm response curve).
    将 LED 光谱与传感器响应重采样到统一细网格后做重叠积分，自动处理采样间隔不一致。

    Parameters
    ----------
    channel : LEDChannel
        The LED. Its measured ``spectrum`` is used if present, otherwise a
        Gaussian line of width ``fwhm_nm`` centered at ``wavelength_nm``.
    response_wavelengths_nm, response_values : np.ndarray
        The sensor spectral response samples (any grid).
    resample_step_nm : float
        Step (nm) of the common grid used for the overlap integral.

    Returns
    -------
    float
        The overlap integral ``∫ I_led(λ) · η(λ) dλ``.
    """
    # Determine the LED's wavelength support.
    if channel.spectrum is not None:
        led_wl = np.asarray(channel.spectrum["wavelengths"], dtype=float)
        led_intensity = np.asarray(channel.spectrum["intensity"], dtype=float)
        led_lo, led_hi = led_wl.min(), led_wl.max()
    else:
        led_wl = led_intensity = None
        # A Gaussian is effectively supported within ~3 FWHM of its center.
        led_lo = channel.wavelength_nm - 3 * channel.fwhm_nm
        led_hi = channel.wavelength_nm + 3 * channel.fwhm_nm

    lo = min(led_lo, response_wavelengths_nm.min())
    hi = max(led_hi, response_wavelengths_nm.max())
    grid = np.arange(lo, hi + resample_step_nm, resample_step_nm)

    # Resample both curves onto the common grid (zero-fill outside support).
    response = np.interp(grid, response_wavelengths_nm, response_values,
                         left=0.0, right=0.0)
    if led_wl is not None:
        emission = np.interp(grid, led_wl, led_intensity, left=0.0, right=0.0)
    else:
        emission = gaussian_emission(grid, channel.wavelength_nm, channel.fwhm_nm)

    return float(np.trapezoid(emission * response, grid))


# ─────────────────────────────────────────────────────────────────────────
# Data structures
# ─────────────────────────────────────────────────────────────────────────
@dataclass
class LEDChannel:
    """One LED / reconstruction channel.

    Attributes
    ----------
    wavelength_nm : float
        Center wavelength (nm); also the key under which the channel's
        reconstructed reflectance is returned.
    frequency_hz : float
        PWM modulation frequency (Hz).
    phase_rad : float
        Initial phase (radians). I-channel = 0, Q-channel = π/2.
    fwhm_nm : float
        Gaussian line width (nm); used only when ``spectrum`` is absent.
    weight : float, optional
        Pre-calibrated channel weight; overrides the spectral computation.
    spectrum : dict, optional
        Measured LED spectrum ``{"wavelengths": [...], "intensity": [...]}``.
    label : str, optional
        Human-readable label.
    """
    wavelength_nm: float
    frequency_hz: float
    phase_rad: float
    fwhm_nm: float = DEFAULT_FWHM_NM
    weight: Optional[float] = None
    spectrum: Optional[dict] = None
    label: Optional[str] = None


@dataclass
class ReconstructionConfig:
    """Full reconstruction configuration (sensor + LEDs + filter)."""
    sample_rate_hz: float
    channels: List[LEDChannel]
    lpf_cutoff_hz: float = DEFAULT_LPF_CUTOFF_HZ
    lpf_order: int = DEFAULT_LPF_ORDER
    integration_time_s: Optional[float] = None
    integration_delay_s: Optional[float] = None
    sensor_response: Optional[dict] = None  # {"wavelengths": [...], "response": [...]}
    resample_step_nm: float = DEFAULT_RESAMPLE_STEP_NM

    @classmethod
    def from_dict(cls, data: dict) -> "ReconstructionConfig":
        """Build a config from a parsed JSON/YAML mapping.

        Expected schema (see ``config/example_config.yaml``)::

            sensor: {fs, integration_time?, integration_delay?, spectral_response?}
            lpf:    {cutoff?, order?}
            leds:   [{wavelength, freq, phase|phase_deg|phase_rad, fwhm?, weight?, spectrum?, label?}, ...]
            resample_step?: float
        """
        sensor = data.get("sensor", {})
        lpf = data.get("lpf", {})

        channels: List[LEDChannel] = []
        for item in data["leds"]:
            channels.append(LEDChannel(
                wavelength_nm=float(item["wavelength"]),
                frequency_hz=float(item["freq"]),
                phase_rad=_parse_phase_rad(item),
                fwhm_nm=float(item.get("fwhm", DEFAULT_FWHM_NM)),
                weight=_optional_float(item.get("weight")),
                spectrum=item.get("spectrum"),
                label=item.get("label"),
            ))

        return cls(
            sample_rate_hz=float(sensor["fs"]),
            channels=channels,
            lpf_cutoff_hz=float(lpf.get("cutoff", DEFAULT_LPF_CUTOFF_HZ)),
            lpf_order=int(lpf.get("order", DEFAULT_LPF_ORDER)),
            integration_time_s=_optional_float(sensor.get("integration_time")),
            integration_delay_s=_optional_float(sensor.get("integration_delay")),
            sensor_response=sensor.get("spectral_response"),
            resample_step_nm=float(data.get("resample_step", DEFAULT_RESAMPLE_STEP_NM)),
        )

    @classmethod
    def from_file(cls, path: str) -> "ReconstructionConfig":
        """Load a config from a ``.json``, ``.yaml`` or ``.yml`` file."""
        with open(path, "r", encoding="utf-8") as handle:
            text = handle.read()
        extension = os.path.splitext(path)[1].lower()
        if extension in (".yaml", ".yml"):
            if not _HAS_YAML:
                raise RuntimeError(
                    "PyYAML is required to read YAML configs "
                    "(`pip install pyyaml`), or use a JSON config.")
            data = yaml.safe_load(text)
        else:
            data = json.loads(text)
        return cls.from_dict(data)


@dataclass
class ReconstructionResult:
    """Output of a reconstruction call.

    Attributes
    ----------
    reflectance : dict
        Maps ``wavelength_nm`` (float) to the reconstructed reflectance series.
    time_s : np.ndarray
        Time axis (s) of the returned series.
    calibration_mode : str
        One of ``"white_reference"``, ``"weights"``, ``"spectral"``, ``"none"``.
    integration_delay_s : float
        Group delay (s) used for the delay-compensated reference.
    meta : dict
        Miscellaneous metadata (sample rate, cutoff, channel count, ...).
    """
    reflectance: Dict[float, np.ndarray]
    time_s: np.ndarray
    calibration_mode: str
    integration_delay_s: float
    meta: dict = field(default_factory=dict)

    @property
    def wavelengths_nm(self) -> List[float]:
        """Sorted list of reconstructed wavelengths (nm)."""
        return sorted(self.reflectance.keys())

    def as_matrix(self) -> Tuple[List[float], np.ndarray]:
        """Return ``(wavelengths_nm, T×C reflectance matrix)``."""
        wavelengths = self.wavelengths_nm
        matrix = np.column_stack([self.reflectance[wl] for wl in wavelengths])
        return wavelengths, matrix


# ─────────────────────────────────────────────────────────────────────────
# Parsing helpers
# ─────────────────────────────────────────────────────────────────────────
def _optional_float(value) -> Optional[float]:
    """Cast to float, preserving ``None``."""
    return None if value is None else float(value)


def _parse_phase_rad(item: dict) -> float:
    """Resolve an LED initial phase to radians.

    Accepts ``phase_rad`` (radians), ``phase_deg`` (degrees), or a bare
    ``phase`` interpreted as **degrees** (the common authoring convention).
    """
    if item.get("phase_rad") is not None:
        return float(item["phase_rad"])
    if item.get("phase_deg") is not None:
        return np.deg2rad(float(item["phase_deg"]))
    if item.get("phase") is not None:
        return np.deg2rad(float(item["phase"]))  # bare 'phase' -> degrees
    return 0.0


# ─────────────────────────────────────────────────────────────────────────
# Reconstructor
# ─────────────────────────────────────────────────────────────────────────
class SpectralReconstructor:
    """IQ phase-encoded lock-in reconstructor.

    Examples
    --------
    >>> rec = SpectralReconstructor.from_config_file("config/example_config.yaml")
    >>> result = rec.reconstruct(signal, white_reference=gray, reference_level=0.5)
    >>> reflectance_850 = result.reflectance[850.0]
    """

    def __init__(self, config: ReconstructionConfig) -> None:
        self.config = config
        # Zero-phase Butterworth low-pass used for every channel's demodulation.
        self._sos = butter(config.lpf_order,
                           config.lpf_cutoff_hz / (config.sample_rate_hz / 2.0),
                           btype="low", output="sos")

    # ----- constructors ----------------------------------------------------
    @classmethod
    def from_config_file(cls, path: str) -> "SpectralReconstructor":
        """Create a reconstructor from a JSON/YAML config file."""
        return cls(ReconstructionConfig.from_file(path))

    @classmethod
    def from_dict(cls, data: dict) -> "SpectralReconstructor":
        """Create a reconstructor from a parsed config mapping."""
        return cls(ReconstructionConfig.from_dict(data))

    # ----- derived quantities ---------------------------------------------
    def group_delay_s(self) -> float:
        """Sensor integration group delay τ (s).

        A boxcar integrator over an exposure window of length ``T`` delays the
        carrier centroid by ``T / 2`` (see the theory doc, boxcar group-delay
        lemma). Resolution order: explicit ``integration_delay_s`` →
        ``integration_time_s / 2`` → full-frame assumption ``(1/fs) / 2``.
        """
        config = self.config
        if config.integration_delay_s is not None:
            return config.integration_delay_s
        exposure_s = (config.integration_time_s
                      if config.integration_time_s is not None
                      else 1.0 / config.sample_rate_hz)
        return 0.5 * exposure_s

    def _reference(self, channel: LEDChannel, time_s: np.ndarray,
                   delay_s: float) -> np.ndarray:
        """Delay-compensated lock-in reference matched to the LED fundamental.

        The unipolar 50%-duty PWM fundamental is ``(2/π) sin(2π f t + φ)``, so
        the matched, delay-compensated reference is ``sin(2π f (t + τ) + φ)``.
        """
        return np.sin(2 * np.pi * channel.frequency_hz * (time_s + delay_s)
                      + channel.phase_rad)

    def _spectral_weights(self) -> Optional[Dict[float, float]]:
        """Per-channel weights from the sensor response, or ``None`` if absent."""
        response = self.config.sensor_response
        if response is None:
            return None
        response_wl = np.asarray(response["wavelengths"], dtype=float)
        response_val = np.asarray(response["response"], dtype=float)
        return {ch.wavelength_nm: spectral_overlap(
                    ch, response_wl, response_val, self.config.resample_step_nm)
                for ch in self.config.channels}

    # ----- main API --------------------------------------------------------
    def reconstruct(self,
                    signal: np.ndarray,
                    time_s: Optional[np.ndarray] = None,
                    white_reference: Optional[np.ndarray] = None,
                    reference_level: float = 1.0,
                    trim_s: float = 0.0) -> ReconstructionResult:
        """Reconstruct per-wavelength reflectance from a captured signal.

        Parameters
        ----------
        signal : np.ndarray
            The captured detector time series (length N).
        time_s : np.ndarray, optional
            Time axis (s). Defaults to ``arange(N) / sample_rate_hz``.
        white_reference : np.ndarray, optional
            A capture of a flat reference target (length N). If given, the
            calibration mode is ``white_reference`` and the output is absolute
            reflectance, independent of LED power or sensor response. The
            target need not be perfectly white — use a gray card and set
            ``reference_level`` to its known reflectance to avoid saturating
            the detector.
        reference_level : float
            Known flat reflectance of the reference target (default 1.0).
        trim_s : float
            Seconds to discard at each end (IIR edge transient).

        Returns
        -------
        ReconstructionResult
        """
        signal = np.asarray(signal, dtype=float).ravel()
        n_samples = signal.size
        sample_rate_hz = self.config.sample_rate_hz

        if time_s is None:
            time_s = np.arange(n_samples) / sample_rate_hz
        else:
            time_s = np.asarray(time_s, dtype=float).ravel()
            if time_s.size != n_samples:
                raise ValueError("len(time_s) must equal len(signal)")

        delay_s = self.group_delay_s()
        mode, white_amplitude, weights = self._resolve_calibration(
            signal, time_s, delay_s, white_reference, n_samples)

        # Demodulate every channel: mix with the matched reference, low-pass,
        # then calibrate to reflectance.
        reflectance: Dict[float, np.ndarray] = {}
        for channel in self.config.channels:
            reference = self._reference(channel, time_s, delay_s)
            low_passed = sosfiltfilt(self._sos, signal * reference)
            if mode == "white_reference":
                reflectance[channel.wavelength_nm] = self._apply_white_reference(
                    low_passed, white_amplitude[channel.wavelength_nm],
                    reference_level, channel.wavelength_nm)
            else:
                # Theorem: LPF[s · ref] = w·R/π  ⇒  R = π·LPF / w.
                reflectance[channel.wavelength_nm] = (
                    np.pi * low_passed / weights[channel.wavelength_nm])

        if trim_s > 0:
            keep = (time_s >= time_s[0] + trim_s) & (time_s <= time_s[-1] - trim_s)
            time_out = time_s[keep]
            reflectance = {wl: series[keep] for wl, series in reflectance.items()}
        else:
            time_out = time_s

        return ReconstructionResult(
            reflectance=reflectance,
            time_s=time_out,
            calibration_mode=mode,
            integration_delay_s=delay_s,
            meta={"sample_rate_hz": sample_rate_hz,
                  "lpf_cutoff_hz": self.config.lpf_cutoff_hz,
                  "lpf_order": self.config.lpf_order,
                  "n_channels": len(self.config.channels)},
        )

    # ----- calibration internals ------------------------------------------
    def _resolve_calibration(self, signal, time_s, delay_s, white_reference,
                             n_samples):
        """Pick the calibration mode and precompute its constants.

        Returns ``(mode, white_amplitude_or_None, weights_or_None)``.
        """
        if white_reference is not None:
            return ("white_reference",
                    self._white_amplitudes(white_reference, time_s, delay_s,
                                           n_samples),
                    None)

        if any(ch.weight is not None for ch in self.config.channels):
            weights = {ch.wavelength_nm: (1.0 if ch.weight is None else ch.weight)
                       for ch in self.config.channels}
            return "weights", None, weights

        weights = self._spectral_weights()
        if weights is not None:
            warnings.warn(
                "Calibrating from the sensor spectral response: output is "
                "RELATIVE reflectance (absolute scale depends on unknown LED "
                "radiant power). Use white_reference for absolute values.")
            return "spectral", None, weights

        warnings.warn(
            "No spectral response, weights, or white reference provided: "
            "returning reflectance up to an unknown per-channel scale.")
        weights = {ch.wavelength_nm: 1.0 for ch in self.config.channels}
        return "none", None, weights

    def _white_amplitudes(self, white_reference, time_s, delay_s,
                          n_samples) -> Dict[float, float]:
        """Per-channel lock-in amplitude of a flat reference capture.

        The reference reflectance is constant, so each channel's lock-in output
        is a single calibration **constant** (≈ w·R_ref/π). A robust interior
        median is used instead of the raw series, which would ring through zero
        at the IIR edges. 取内区中位数而非逐点序列，避免边缘振铃过零。
        """
        white = np.asarray(white_reference, dtype=float).ravel()
        if white.size != n_samples:
            raise ValueError("white_reference must match signal length")
        lo, hi = int(0.1 * n_samples), int(0.9 * n_samples)
        amplitudes: Dict[float, float] = {}
        for channel in self.config.channels:
            reference = self._reference(channel, time_s, delay_s)
            series = sosfiltfilt(self._sos, white * reference)
            amplitudes[channel.wavelength_nm] = float(np.median(series[lo:hi]))
        return amplitudes

    @staticmethod
    def _apply_white_reference(low_passed, white_amplitude, reference_level,
                               wavelength_nm) -> np.ndarray:
        """Divide a channel's lock-in output by its white-reference constant."""
        if abs(white_amplitude) < 1e-12:
            warnings.warn(f"white reference ~0 for {wavelength_nm} nm; "
                          "channel left uncalibrated")
            return np.full_like(low_passed, np.nan)
        # R = R_ref · LPF[s·ref] / LPF[s_ref·ref]; w and π cancel exactly.
        return reference_level * low_passed / white_amplitude


# ─────────────────────────────────────────────────────────────────────────
# CLI: python spectral_reconstruction.py <config> <signal> [white_reference]
# ─────────────────────────────────────────────────────────────────────────
def _load_series(path: str) -> np.ndarray:
    """Load a 1-D series from ``.npy`` or a delimited text file."""
    if path.lower().endswith(".npy"):
        return np.load(path)
    return np.loadtxt(path, delimiter=",")


def _main() -> None:
    import sys
    if len(sys.argv) < 3:
        print("usage: python spectral_reconstruction.py <config> "
              "<signal.npy|csv> [white_reference.npy|csv]")
        raise SystemExit(1)

    reconstructor = SpectralReconstructor.from_config_file(sys.argv[1])
    signal = _load_series(sys.argv[2])
    white = _load_series(sys.argv[3]) if len(sys.argv) > 3 else None
    result = reconstructor.reconstruct(signal, white_reference=white, trim_s=1.0)

    print(f"calibration mode : {result.calibration_mode}")
    print(f"integration delay: {result.integration_delay_s * 1e3:.3f} ms")
    print(f"channels         : {result.meta['n_channels']}")
    for wavelength_nm in result.wavelengths_nm:
        series = result.reflectance[wavelength_nm]
        print(f"  {wavelength_nm:>7.1f} nm  mean={np.nanmean(series):.4f}  "
              f"min={np.nanmin(series):.4f}  max={np.nanmax(series):.4f}")


if __name__ == "__main__":
    _main()
