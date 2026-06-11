"""
spectral_reconstruction.py
==========================
Config-driven multispectral reconstruction via IQ phase-encoded lock-in
demodulation.

配置驱动的多光谱重建接口（IQ 相位编码锁相解调）。

This module turns the algorithm validated in ``iq_sensing_system.py`` into a
reusable interface:

    1. Load a configuration file (JSON or YAML) describing the **sensor** and
       the **LEDs** (frequency, phase, wavelength, optional spectrum / weight).
    2. Feed it the captured 1-D time-series signal from the detector.
    3. Get back the reconstructed time-varying reflectance for every LED
       wavelength.

Design notes
------------
* **The sensor spectral response is OPTIONAL.** In practice it is usually
  unknown. Three calibration modes are supported, in order of practicality:

      - ``white_reference``  (recommended) — supply one capture of a flat white
        target; absolute reflectance is obtained with no knowledge of LED power
        or sensor response.
      - ``weights``          — per-channel weights known from calibration.
      - ``spectral``         — compute weights from the (optional) sensor
        spectral response + LED spectra. Gives *relative* reflectance only
        (absolute scale depends on unknown LED radiant power).
      - ``none``             — no calibration; returns reflectance up to an
        unknown per-channel scale (still correct in shape / dynamics).

* **Mismatched wavelength sampling is handled.** A datasheet LED spectrum may be
  sampled every 5 nm while the sensor response is sampled every 1 nm (or 10 nm).
  Both are interpolated onto a common fine grid (``resample_step`` nm, default
  1 nm) with zero-fill outside their support before the overlap integral.

Dependencies: numpy, scipy. YAML configs additionally require PyYAML.
依赖：numpy、scipy；使用 YAML 配置时还需 PyYAML。
"""

from __future__ import annotations

import json
import os
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.signal import butter, sosfiltfilt

try:
    import yaml  # optional, only needed for .yaml/.yml configs
    _HAS_YAML = True
except Exception:  # pragma: no cover
    _HAS_YAML = False


# ─────────────────────────────────────────────────────────────────────────
# Data structures
# ─────────────────────────────────────────────────────────────────────────
@dataclass
class LEDChannel:
    """One LED / reconstruction channel.

    wavelength : center wavelength (nm), used as the channel key
    freq       : modulation frequency (Hz)
    phase      : initial phase (radians); I-channel=0, Q-channel=pi/2
    fwhm       : Gaussian line-width (nm), used only if ``spectrum`` is absent
    weight     : optional pre-calibrated channel weight (overrides spectral calc)
    spectrum   : optional measured LED spectrum {'wavelengths': [...], 'intensity': [...]}
    label      : optional human label
    """
    wavelength: float
    freq: float
    phase: float
    fwhm: float = 20.0
    weight: Optional[float] = None
    spectrum: Optional[dict] = None
    label: Optional[str] = None


@dataclass
class ReconConfig:
    """Full reconstruction configuration."""
    fs: float                                   # sensor sample rate (Hz)
    leds: list                                  # list[LEDChannel]
    cutoff: float = 3.5                          # lock-in LPF cutoff (Hz)
    order: int = 4                               # Butterworth order
    integration_time: Optional[float] = None     # sensor exposure / integration (s)
    integration_delay: Optional[float] = None     # explicit boxcar group delay (s); overrides
    sensor_response: Optional[dict] = None        # optional {'wavelengths':[...], 'response':[...]}
    resample_step: float = 1.0                    # common-grid step (nm) for spectral overlap

    # ----- construction helpers -------------------------------------------
    @classmethod
    def from_dict(cls, d: dict) -> "ReconConfig":
        sensor = d.get("sensor", {})
        lpf = d.get("lpf", {})

        leds = []
        for i, item in enumerate(d["leds"]):
            phase = _parse_phase(item)
            leds.append(LEDChannel(
                wavelength=float(item["wavelength"]),
                freq=float(item["freq"]),
                phase=phase,
                fwhm=float(item.get("fwhm", 20.0)),
                weight=(None if item.get("weight") is None else float(item["weight"])),
                spectrum=item.get("spectrum"),
                label=item.get("label"),
            ))

        return cls(
            fs=float(sensor["fs"]),
            leds=leds,
            cutoff=float(lpf.get("cutoff", 3.5)),
            order=int(lpf.get("order", 4)),
            integration_time=_opt_float(sensor.get("integration_time")),
            integration_delay=_opt_float(sensor.get("integration_delay")),
            sensor_response=sensor.get("spectral_response"),
            resample_step=float(d.get("resample_step", 1.0)),
        )

    @classmethod
    def from_file(cls, path: str) -> "ReconConfig":
        with open(path, "r", encoding="utf-8") as fh:
            text = fh.read()
        ext = os.path.splitext(path)[1].lower()
        if ext in (".yaml", ".yml"):
            if not _HAS_YAML:
                raise RuntimeError(
                    "PyYAML is required to read YAML configs. "
                    "Install it (`pip install pyyaml`) or use a JSON config."
                )
            d = yaml.safe_load(text)
        else:
            d = json.loads(text)
        return cls.from_dict(d)


@dataclass
class ReconResult:
    """Result of a reconstruction call."""
    reflectance: dict                # {wavelength(float): np.ndarray}
    t: np.ndarray                    # time axis (s) of the returned series
    calibration_mode: str            # 'white_reference' | 'weights' | 'spectral' | 'none'
    integration_delay: float         # group delay used (s)
    meta: dict = field(default_factory=dict)

    @property
    def wavelengths(self):
        return sorted(self.reflectance.keys())

    def as_matrix(self):
        """Return (wavelengths, T x C reflectance matrix)."""
        wls = self.wavelengths
        M = np.column_stack([self.reflectance[w] for w in wls])
        return wls, M


# ─────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────
def _opt_float(x):
    return None if x is None else float(x)


def _parse_phase(item: dict) -> float:
    """Accept phase as radians ('phase_rad'), degrees ('phase_deg'),
    or bare 'phase' interpreted as DEGREES (the common convention)."""
    if "phase_rad" in item and item["phase_rad"] is not None:
        return float(item["phase_rad"])
    if "phase_deg" in item and item["phase_deg"] is not None:
        return np.deg2rad(float(item["phase_deg"]))
    if "phase" in item and item["phase"] is not None:
        return np.deg2rad(float(item["phase"]))   # bare 'phase' -> degrees
    return 0.0


def gaussian_line(grid: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    """Normalized-peak Gaussian emission line on ``grid`` (nm)."""
    return np.exp(-4.0 * np.log(2.0) * ((grid - center) / fwhm) ** 2)


def overlap_integral(led: LEDChannel,
                     resp_wl: np.ndarray,
                     resp_val: np.ndarray,
                     step: float = 1.0) -> float:
    """Spectral overlap of one LED with the sensor response.

    Resamples BOTH curves onto a common fine grid (handles different native
    wavelength sampling intervals) and integrates their product.

    将 LED 光谱与传感器响应重采样到统一细网格后做重叠积分，
    自动处理两者波长采样间隔不一致的情况。
    """
    # determine LED support
    if led.spectrum is not None:
        lw = np.asarray(led.spectrum["wavelengths"], float)
        li = np.asarray(led.spectrum["intensity"], float)
        led_lo, led_hi = lw.min(), lw.max()
    else:
        # Gaussian effectively supported within ~3*FWHM of the center
        lw = li = None
        led_lo, led_hi = led.wavelength - 3 * led.fwhm, led.wavelength + 3 * led.fwhm

    lo = min(led_lo, resp_wl.min())
    hi = max(led_hi, resp_wl.max())
    grid = np.arange(lo, hi + step, step)

    R = np.interp(grid, resp_wl, resp_val, left=0.0, right=0.0)
    if lw is not None:
        L = np.interp(grid, lw, li, left=0.0, right=0.0)
    else:
        L = gaussian_line(grid, led.wavelength, led.fwhm)

    return float(np.trapezoid(L * R, grid))


# ─────────────────────────────────────────────────────────────────────────
# Reconstructor
# ─────────────────────────────────────────────────────────────────────────
class SpectralReconstructor:
    """IQ phase-encoded lock-in reconstructor.

    Examples
    --------
    >>> rec = SpectralReconstructor.from_config_file("config/example_config.yaml")
    >>> result = rec.reconstruct(signal)                 # uncalibrated / spectral
    >>> result = rec.reconstruct(signal, white_reference=white_signal)  # absolute
    >>> refl_850 = result.reflectance[850.0]
    """

    def __init__(self, config: ReconConfig):
        self.cfg = config
        self._sos = butter(config.order,
                           config.cutoff / (config.fs / 2.0),
                           btype="low", output="sos")

    # ----- constructors ----------------------------------------------------
    @classmethod
    def from_config_file(cls, path: str) -> "SpectralReconstructor":
        return cls(ReconConfig.from_file(path))

    @classmethod
    def from_dict(cls, d: dict) -> "SpectralReconstructor":
        return cls(ReconConfig.from_dict(d))

    # ----- derived quantities ---------------------------------------------
    def group_delay(self) -> float:
        """Sensor integration group delay tau (s).

        Priority: explicit ``integration_delay`` > ``integration_time``/2 >
        full-frame assumption (1/fs)/2.

        A boxcar integrator over an exposure window of length T has its
        centroid delayed by T/2 (see the theory document, Lemma on boxcar
        group delay).
        """
        cfg = self.cfg
        if cfg.integration_delay is not None:
            return cfg.integration_delay
        t_int = cfg.integration_time if cfg.integration_time is not None else 1.0 / cfg.fs
        return 0.5 * t_int

    def _spectral_weights(self) -> Optional[dict]:
        """Per-channel weights from sensor response, or None if unavailable."""
        resp = self.cfg.sensor_response
        if resp is None:
            return None
        rw = np.asarray(resp["wavelengths"], float)
        rr = np.asarray(resp["response"], float)
        return {led.wavelength: overlap_integral(led, rw, rr, self.cfg.resample_step)
                for led in self.cfg.leds}

    def _reference(self, led: LEDChannel, t: np.ndarray, tau: float) -> np.ndarray:
        """Delay-compensated lock-in reference matched to the LED fundamental.

        The unipolar 50%-duty PWM fundamental is (2/pi) sin(2*pi f t + phi),
        so the matched reference is sin(2*pi f (t + tau) + phi).
        """
        return np.sin(2 * np.pi * led.freq * (t + tau) + led.phase)

    # ----- main API --------------------------------------------------------
    def reconstruct(self,
                    signal: np.ndarray,
                    t: Optional[np.ndarray] = None,
                    white_reference: Optional[np.ndarray] = None,
                    reference_level: float = 1.0,
                    trim: float = 0.0) -> ReconResult:
        """Reconstruct per-wavelength reflectance from the captured signal.

        Parameters
        ----------
        signal : 1-D array, the captured detector time series (length N).
        t      : optional time axis (s). Defaults to ``arange(N)/fs``.
        white_reference : optional capture of a flat reference target (length N).
                 If given, calibration mode is ``white_reference`` and the
                 output is absolute reflectance regardless of LED power or
                 sensor response. The reference need not be perfectly white --
                 use a gray target and set ``reference_level`` to its known
                 reflectance to avoid saturating the detector.
        reference_level : known flat reflectance of the reference target
                 (default 1.0 = ideal white).
        trim   : seconds to discard at each end (IIR edge transient).

        Returns
        -------
        ReconResult
        """
        signal = np.asarray(signal, float).ravel()
        N = signal.size
        fs = self.cfg.fs
        if t is None:
            t = np.arange(N) / fs
        else:
            t = np.asarray(t, float).ravel()
            if t.size != N:
                raise ValueError("len(t) must equal len(signal)")

        tau = self.group_delay()

        # ---- decide calibration mode ------------------------------------
        white_amp = None
        weights = None
        if white_reference is not None:
            mode = "white_reference"
            wref = np.asarray(white_reference, float).ravel()
            if wref.size != N:
                raise ValueError("white_reference must match signal length")
            # The white target has constant reflectance (R=1), so its lock-in
            # output is a single calibration CONSTANT per channel (= w_k/pi).
            # Take a robust scalar (interior median) rather than the raw series,
            # which would ring through zero at the IIR edges.
            white_amp = {}
            lo, hi = int(0.1 * N), int(0.9 * N)
            for led in self.cfg.leds:
                ref = self._reference(led, t, tau)
                wa = sosfiltfilt(self._sos, wref * ref)
                white_amp[led.wavelength] = float(np.median(wa[lo:hi]))
        elif any(led.weight is not None for led in self.cfg.leds):
            mode = "weights"
            weights = {led.wavelength: (1.0 if led.weight is None else led.weight)
                       for led in self.cfg.leds}
        else:
            weights = self._spectral_weights()
            if weights is not None:
                mode = "spectral"
                warnings.warn(
                    "Calibrating from sensor spectral response: output is "
                    "RELATIVE reflectance (absolute scale depends on unknown "
                    "LED radiant power). Use white_reference for absolute values."
                )
            else:
                mode = "none"
                weights = {led.wavelength: 1.0 for led in self.cfg.leds}
                warnings.warn(
                    "No spectral response, weights, or white reference provided: "
                    "returning reflectance up to an unknown per-channel scale."
                )

        # ---- demodulate every channel -----------------------------------
        refl = {}
        for led in self.cfg.leds:
            ref = self._reference(led, t, tau)
            lp = sosfiltfilt(self._sos, signal * ref)
            if mode == "white_reference":
                wa = white_amp[led.wavelength]   # scalar calibration constant
                if abs(wa) < 1e-12:
                    warnings.warn(f"white reference ~0 for {led.wavelength} nm; "
                                  "channel left uncalibrated")
                    refl[led.wavelength] = np.full_like(lp, np.nan)
                else:
                    refl[led.wavelength] = reference_level * lp / wa
            else:
                refl[led.wavelength] = np.pi * lp / weights[led.wavelength]

        # ---- trim edges --------------------------------------------------
        if trim > 0:
            keep = (t >= t[0] + trim) & (t <= t[-1] - trim)
            t_out = t[keep]
            refl = {wl: r[keep] for wl, r in refl.items()}
        else:
            t_out = t

        return ReconResult(
            reflectance=refl,
            t=t_out,
            calibration_mode=mode,
            integration_delay=tau,
            meta={"fs": fs, "cutoff": self.cfg.cutoff, "order": self.cfg.order,
                  "n_channels": len(self.cfg.leds)},
        )


# ─────────────────────────────────────────────────────────────────────────
# CLI: python spectral_reconstruction.py <config> <signal.npy|csv> [white.npy]
# ─────────────────────────────────────────────────────────────────────────
def _load_series(path: str) -> np.ndarray:
    if path.lower().endswith(".npy"):
        return np.load(path)
    return np.loadtxt(path, delimiter=",")


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("usage: python spectral_reconstruction.py <config> "
              "<signal.npy|csv> [white_reference.npy|csv]")
        raise SystemExit(1)

    rec = SpectralReconstructor.from_config_file(sys.argv[1])
    sig = _load_series(sys.argv[2])
    white = _load_series(sys.argv[3]) if len(sys.argv) > 3 else None
    res = rec.reconstruct(sig, white_reference=white, trim=1.0)

    print(f"calibration mode : {res.calibration_mode}")
    print(f"integration delay: {res.integration_delay*1e3:.3f} ms")
    print(f"channels         : {res.meta['n_channels']}")
    for wl in res.wavelengths:
        r = res.reflectance[wl]
        print(f"  {wl:>7.1f} nm  mean={np.nanmean(r):.4f}  "
              f"min={np.nanmin(r):.4f}  max={np.nanmax(r):.4f}")
