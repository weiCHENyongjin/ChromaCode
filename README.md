# ChromaCode

**Multi-spectral video through a single detector.**

ChromaCode recovers reflectance across many spectral bands from a single detector. Instead of a filter wheel, a camera array, or a prism, a scene is captured through one sensor and the spectral bands are recovered from it in software. The same recovery runs per pixel, so an ordinary area sensor yields multi-spectral images and video.

![Reconstruction overview](assets/overview.png)

*Ground truth, reconstruction, and per-pixel error for the same scene. Left to right: reference reflectance, recovered reflectance, absolute error.*

## What it does

- **One detector, many bands.** Multiple spectral channels share a single acquisition path; the spectrometer is software, not optics.
- **Runs per pixel.** Recovery runs across the whole frame, so a full-resolution sensor produces multi-spectral video rather than a single-point reading.
- **Batch and live output.** Analysis can run over a recording for maximum accuracy, or continuously for live output.
- **Sized for real hardware.** Reconstruction stays within interactive time budgets at full HD.

## Results

Reconstruction quality on a slowly varying scene. The error map stays near zero across the frame, and per-band error is uniform rather than concentrated in a few channels.

![Per-band reflectance error](assets/per-band-error.png)

*Reconstruction RMSE for each spectral band. Every band stays well below the 0.05 reference line.*

Spectra are recovered, not merely classified. Extracting individual pixels and comparing recovered curves against ground truth shows the shape of each spectrum is preserved.

![Spectral curve comparison](assets/spectral-curves.png)

*Recovered vs. reference spectra sampled at three positions across the scene.*

## Reconstruction in motion

A moving scene, reconstructed band by band. The trajectory is drawn to make motion visible; the error panel tracks the reconstruction against ground truth frame by frame.

![Dynamic reconstruction](assets/reconstruction.gif)

## Streaming output

The streaming path updates every sample rather than waiting for a full batch. The left panel shows the mixed acquisition signal at a single pixel; the right panel shows each band converging to its reflectance value.

![Streaming reconstruction](assets/streaming-demo.gif)

## Applications

Spectral video is useful wherever spectral information matters but a fixed filter stack or camera array is too heavy, too slow, or too expensive.

![Application areas](assets/applications.png)

Agricultural sensing, industrial inspection and material sorting, medical and life-science imaging, cultural-heritage analysis, low-cost robotic spectral vision, and environmental monitoring.

## Status

ChromaCode is under active development. Source code and technical documentation are not published in this repository; this repository hosts the project overview, visual results, and demos only.

## Contact

For collaboration or licensing enquiries, open an issue or reach out directly.
