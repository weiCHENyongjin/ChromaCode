# Code conventions / 代码规范

This project follows a single, strict convention so the code stays easy to read
for newcomers. 本项目遵循统一、严格的代码规范，便于他人阅读理解。

## Language / 语言

- **Code identifiers and docstrings: English.** A single language keeps the API
  consistent and internationally readable.
- **Inline comments: English, with a short Chinese gloss where a concept is
  subtle.** 关键概念可加简短中文说明。
- User-facing prose (READMEs) is bilingual; see `README.md` / `README.zh-CN.md`.

## Style / 风格

- **PEP 8** layout; **PEP 257** docstrings in **NumPy style** (`Parameters`,
  `Returns`, `Raises`). Lines ≤ 88 chars where practical.
- **Type hints** on every public function, method, and dataclass field.
- Naming:
  - `snake_case` — functions, methods, variables, module-level config keys.
  - `PascalCase` — classes and dataclasses.
  - `UPPER_SNAKE_CASE` — module-level constants.
- **Units in names.** Physical quantities carry their unit as a suffix:
  `frequency_hz`, `wavelength_nm`, `integration_time_s`, `phase_rad`. Domain
  symbols that match the theory document (`t`, `fs`, `f`) are kept short *inside
  local numeric code* but documented.
- **No magic numbers** in logic — promote to a named constant or config field
  with a comment explaining its origin.
- **Single source of truth.** The lock-in demodulation lives only in
  `spectral_reconstruction.py`; the simulator reuses it rather than duplicating.

## Module layout / 模块结构

Each module is ordered: module docstring → imports → constants → data structures
(dataclasses) → pure helper functions → classes → `main()` / CLI guard.

## Testing / 测试

`tests/test_reconstruction.py` drives the public API end-to-end and asserts the
validated accuracy (mean RMSE ≈ 0.044). Run it after any change:

```bash
python tests/test_reconstruction.py
```
