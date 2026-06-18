# Refactor design — re-base `dssatcalibrator` on the shared stack

**Status:** implemented behind feature flags (2026-06-18). **Drafted:** 2026-06-18.
**Goal:** make `dssatcalibrator` consume the workspace's shared layers — `dssatutils`
(weather/soil acquisition), `dssatengine` (DSSAT run/parse), the central `DSSAT48`
install (binary + genotype + weather + soil), and the shared `dssat_templates/` FileX
templates — **to the extent the calibration use case allows**, instead of carrying its
own forks. Aligns the repo with `AGENTS.md` principles 1–2 (applications consume shared
layers; don't fork engine/utils logic).

> **Headline finding.** A *full* re-basing is **not** possible, because calibration needs
> three things the shared stack does not currently serve: (a) **daily PlantGro** time-series
> and (b) **Evaluate.OUT** simulated-vs-measured pairs (the engine emits only end-of-season
> `summary.csv`), and (c) **per-run genotype-coefficient perturbation** (the engine assembles
> a FileX from coordinates; it has no notion of editing `.CUL/.ECO/.SPE`). So the realistic
> target is a **partial re-basing**: delegate *execution, acquisition, templates, and the
> DSSAT48 layout*; **retain** the calibration-specific *perturbation writers, output parsers,
> objective, and engines*.

---

## 1. Where `dssatcalibrator` forks the shared stack today

| Axis | Current (self-contained) | File | Shared equivalent |
|---|---|---|---|
| DSSAT execution | `subprocess.run([exe, model, "B", batch], stdout=DEVNULL, check=False)` | `spawn.py:196` | `dssatengine.run_dssat` (logged, raises on non-zero) |
| Batch file | `write_dssbatch` | `spawn.py:52` | `dssatengine.write_dssbatch` |
| Treatments | `parse_treatments` | `spawn.py:36` | `dssatengine.normalize_treatment_list` |
| Output parse | `parse_plantgro` / `parse_evaluate` / `parse_summary` | `dssat_io.py` | **none** — engine reads `summary.csv` only |
| Weather acquisition | standalone NASA POWER (`weather.py`) | `weather.py` | `dssatutils.process_weather_*` |
| Soil acquisition | **none** (uses install `.SOL`) | — | `dssatutils.process_soils_*` |
| Weather gap-fill / horizon | `fill_gap` | `weather.py` | `dssatengine.extend_weather_repeat_single_ignore_partial` (partial overlap) |
| Genotype files | read `DSSAT48/Genotype/<stem>.{CUL,ECO,SPE}` | `spawn.py:99` | **already aligned** ✅ |
| Soil/weather resolution | DSSAT resolves from install via `DSSATPRO.V48` | `spawn.py` | **already aligned** ✅ |
| Experiment FileX | read real FileX from `source.hemp_dir` | `spawn.py:124` | `dssat_templates/` (for *synthesised* FileX only) |
| Coefficient/mgmt perturbation | `edit_cultivar/ecotype/species/filex/soil/weather` | `writers.py` | **none** — calibration-unique |

`DEPENDENCIES.md` previously recorded the calibrator as *"does not use `dssatutils`/`dssatengine`
— own config-driven DSSAT wrapper."* This refactor changes that line.

---

## 2. What the shared layers actually expose (verified)

### 2.1 `dssatengine` (`python/dssatengine/engine.py`)
- **Run unit is grid-point + summary-CSV oriented.** `_run_simulation(ID, points_row, …,
  template_file_name, run_mode, …)` builds a FileX from a **template + coordinates**, runs
  DSSAT in mode `A`/`Q`, and reads **`summary.csv`** (requires the FileX OUTPUTS line to end
  in `FMOPT='C'`). It assembles a fixed end-of-season result frame (`top_weight_kg_ha`,
  `final_grain_kg_ha`, soil C/N, …). **It does not read `PlantGro` (daily) or `Evaluate.OUT`.**
- **Reusable, low-level helpers** (public as of `dssatengine@v0.3.0`):
  - `run_dssat(run_dir, exe, run_mode_flag, filex, model=None, timeout=None)` — robust
    executor: resolves `exe`, captures stdout/stderr to a log, **raises `RuntimeError` on
    non-zero exit** (strictly better than the calibrator's old silent `DEVNULL`/`check=False`).
  - `write_dssbatch`, `write_dssbatch_sequence`, `normalize_treatment_list`.
  - `extend_weather_repeat_single_ignore_partial` — weather extension (relevant to the
    forecast horizon).
- Private Python aliases remain for older consumers that imported the pre-`v0.3.0` names.

### 2.2 `dssatutils` (`python/dssatutils/__init__.py`)
- Public `process_weather_*` / `process_soils_*` (NASA POWER, AgERA5, Daymet, GridMET,
  Open-Meteo, ERA5-Land, …; SSURGO, SoilGrids, POLARIS, gNATSGO, …). Each is **batch /
  whole-year / geopandas-`GeoDataFrame`-driven** and writes `.WTH`/`.SOL` to an `output_dir`.
- `process_weather_nasapower(shapefile, start_year, end_year, output_dir, id_col, lat_col,
  lon_col, n_cores, log_file)` already handles current-year partial-data (caps at today−2).
- **Mismatch for in-season:** no single-point convenience, no sub-year window, no
  `climatology`/`persistence` gap-fill or `horizon` extension — those stay in the calibrator.
- **Dependency weight:** importing a `dssatutils` submodule pulls geopandas/xarray/rasterio.
  Must be an **optional extra**, not a hard dep, to keep `pip install dssatcalibrator` light.

### 2.3 `DSSAT48` install + `dssat_templates/`
- `DSSAT48/` holds the **binary**, **`Genotype/` (`.CUL/.ECO/.SPE`)**, **`Weather/`**,
  **`Soil/`**. The calibrator already reads genotype/soil/weather from `dssat_dir` ✅.
- **Experiment FileX templates do *not* live in `DSSAT48`** — they live in the shared
  **`DSSAT_Gridded_Run_Tutorial/dssat_templates/`**, resolved via `DSSAT_TEMPLATE_DIR` /
  `template_dir` (the gridded pipeline: `dssat_main_pipeline.py:284-291`). The calibrator
  instead reads **real** experiment FileX (with treatments + observations) from
  `source.hemp_dir`. Both are legitimate — templates are for *synthesising* a FileX from
  coordinates; calibration against real trials needs the real FileX.

---

## 3. Target architecture (partial re-basing)

```
                         dssatcalibrator  (calibration science — STAYS)
   samplers · spaces · priors · objective · engines/ · forecast · diagnostics · scaffold
        │ proposes θ                                   ▲ daily PlantGro + Evaluate pairs
        ▼                                              │ (parsed by dssat_io — STAYS)
   writers.py  (perturb .CUL/.ECO/.SPE/FileX — STAYS, calibration-unique)
        │
        ▼ run dir materialised
   ┌─────────────────────────── DELEGATE ───────────────────────────┐
   │ dssatengine.run_dssat (+ write_dssbatch, normalize_treatments)  │  ← execution
   │ dssatutils.process_weather_* / process_soils_*                  │  ← acquisition (new sites)
   │ dssat_templates/  +  engine FileX assembly                      │  ← synthesised experiments
   │ DSSAT48/  (binary · Genotype · Weather · Soil)                  │  ← single source of inputs
   └─────────────────────────────────────────────────────────────────┘
```

Two experiment-setup (“Level A”) paths, selected per `exp_id`:
- **Real experiment** (has a shipped FileX + FileA/FileT): read the real FileX from a
  configured experiments dir — *unchanged*. This is the primary calibration path.
- **Synthesised experiment** (new site/cultivar/species from coordinates + planting date):
  use `dssat_templates/` + the engine's template→FileX assembly + `dssatutils` weather/soil.
  This is where the new-crop and bare-site use cases finally consume the shared Level-A.

The per-spawn (“Level B”) flow keeps the calibrator's perturbation + parsing, but swaps the
**executor** for the engine's.

---

## 4. Refactor plan — phased, with upstream prerequisites

### Phase 0 — upstream prerequisites (foundation-first, per AGENTS.md §2)
- **`dssatengine`:** promote `_run_dssat` → `run_dssat`, `_write_dssbatch` → `write_dssbatch`,
  `_normalize_treatment_list` → `normalize_treatment_list`, and re-export
  `extend_weather_repeat_single_ignore_partial`. Mirror in R for parity. Tag `dssatengine@v0.3.0`.
- **(Optional) `dssatengine` parsers:** add `parse_plantgro`/`parse_evaluate` so the daily +
  sim-vs-meas readers become shared too. If declined, the calibrator keeps `dssat_io.py`
  (documented divergence — the engine's CSV-summary contract genuinely doesn't cover it).
- **`dssatcalibrator/pyproject.toml`:** add pinned, **optional** extras
  `[shared]` → `dssatengine @ git+…@v0.3.0`; `[acquire]` → `dssatutils @ …@v0.2.0` (+ geopandas).
  Never `@main`/editable in a committed manifest (principle 3).

### Phase 1 — formalise the `DSSAT48` layer (low risk, no upstream dep)
- Single config block for the install: `dssat_dir` → binary, `Genotype/`, `Weather/`, `Soil/`.
- Replace any implicit assumptions in `spawn.py` with explicit `dssat_dir`-rooted paths
  (already 90% there). Document `DSSAT48` as the one source of binary + genotype + resolved
  weather/soil. **Behaviour-preserving.**

### Phase 2 — execute via `dssatengine` (the highest-value reuse)
- Add `execution.backend: native | dssatengine` (default `native` initially).
- When `dssatengine`: in `spawn.spawn_and_run`, after writers run, call
  `dssatengine.run_dssat(run_dir, exe, "B", model=crop["model"])` and `write_dssbatch` instead of the local
  subprocess + `write_dssbatch`. **Keep** `dssat_io.parse_plantgro/parse_evaluate` for outputs.
- Gain: logged, fail-loud execution (no more silent `DEVNULL`); shared treatment handling.
- **Parity test** (slow tier): native vs `dssatengine` backend must produce identical PlantGro
  and identical objective score on the hemp example. Flip the default to `dssatengine` once green.

### Phase 3 — acquire weather/soil via `dssatutils` (new-site path only)
- New `weather.provider: dssatutils` that builds a 1-row `GeoDataFrame` and calls
  `process_weather_nasapower(...)` (or any `process_weather_*`), writing `.WTH` into the run/cache.
- Add a soil acquisition path (`process_soils_ssurgo/soilgrids`) for sites with no shipped `.SOL`.
- **Keep** the in-season layer (`fill_gap` horizon/`climatology`, `obs_operator`, anchor) on top —
  optionally back `fill_gap` with `engine.extend_weather_repeat_single_ignore_partial`.
- The standalone NASA POWER fetch in `weather.py` becomes the **zero-dependency fallback**
  (`provider: nasa_power`); `provider: dssatutils` becomes the recommended path when the
  `[acquire]` extra is installed. Removes the documented fork.

### Phase 4 — align templates with `dssat_templates/`
- For **synthesised** experiments, read FileX templates from `DSSAT_TEMPLATE_DIR` /
  `template_dir` (default the shared `DSSAT_Gridded_Run_Tutorial/dssat_templates/`), and use
  the engine's template→FileX assembly + the calibrator's placeholder injection.
- Point `scaffold.py`'s cloned `.CUL/.ECO/.SPE` output at the same template-dir convention so
  a scaffolded new crop is discoverable by both the calibrator and the gridded engine.
- **Real-experiment calibration is unchanged** (it needs the real FileX, not a template).

---

## 5. What deliberately stays in `dssatcalibrator` (and why)

| Component | Why it can't/shouldn't be delegated |
|---|---|
| `writers.py` (perturb `.CUL/.ECO/.SPE/FileX/.SOL/.WTH`) | Per-run coefficient perturbation is the calibrator's core; absent from the engine. |
| `dssat_io.py` (`parse_plantgro`, `parse_evaluate`) | Engine emits only end-of-season `summary.csv`; calibration needs **daily** series + **sim-vs-meas** pairs. (Unless Phase 0 adds them upstream.) |
| `objective.py`, `engines/`, `forecast.py`, `diagnostics.py`, samplers/spaces/priors | Pure calibration science — not in scope for any shared layer. |
| In-season gap-fill/horizon/anchor/obs-operator | Nowcast-specific; no shared equivalent. |

This boundary is the honest answer to "to the extent possible": **reuse the plumbing
(execute, acquire, assemble, locate inputs); keep the calibration brain.**

---

## 6. The output-mode decision (key design fork)

The engine requires **`FMOPT='C'` (CSV: `summary.csv`)**; the calibrator parses classic
**`.OUT`** files (`PlantGro.OUT`, `Evaluate.OUT`). Two options:

- **(A — recommended) Reuse only the engine's *executor*, keep `.OUT` parsing.** Minimal
  coupling; the calibrator keeps full control of daily + sim-vs-meas extraction; no dependence
  on the engine's evolving summary schema. Cost: the parsers remain a (documented) divergence.
- **(B) Switch the calibrator to CSV mode and extend the engine to parse `plantgro.csv` /
  `evaluate.csv`.** Maximises sharing but couples the calibrator to the engine's output
  contract and requires upstream parser work + verifying CSV mode still emits Evaluate pairs.

Recommend **A** now, leaving **B** as a later consolidation if the engine grows daily parsers.

---

## 7. Dependency, parity & governance impact

- **Pins:** add `dssatutils@v0.2.0` / `dssatengine@v0.3.0` (post-Phase-0 tag) as **optional
  extras**; refresh any lockfile; update `DEPENDENCIES.md` (calibrator row: *"does not use"* →
  *"consumes `dssatengine` for execution; `dssatutils` for acquisition (optional extra)"*).
- **R/Python parity (principle 5):** the calibrator is Python-only today; it stays Python-only
  (documented). Upstream engine API promotions in Phase 0 must be mirrored R↔Python.
- **Caches (principle 6):** `weather_cache/` already git-ignored; engine run logs land under
  `results/` (ignored).
- **Portable paths (principle 7):** all install/template locations via config/env
  (`dssat_dir`, `DSSAT_TEMPLATE_DIR`) — no personal absolute paths.
- **Docs:** update `CONCEPT.md` §2 (position in stack) + §0 table, `README.md`, and the
  `weather.py` divergence note (Phase 3 removes the fork).

---

## 8. Risk & migration safety

- **Feature-flagged** (`execution.backend`, `weather.provider`) with the **native path retained
  as fallback** → instant rollback, no flag-day.
- **Parity tests** (native vs shared backend) gate each phase; behaviour-preserving by design.
- **Heavy deps isolated** behind extras → light default install unaffected.
- **Upstream-first** ordering avoids patching a consumer copy of engine logic.

---

## 9. Recommendation (priority order)

1. **Phase 1 (DSSAT48 formalisation)** — do unconditionally; trivial, no upstream dep.
2. **Phase 0 + Phase 2 (execute via `dssatengine`)** — **highest value**: removes the executor
   fork, gains fail-loud logging, needs only a small public-API promotion upstream.
3. **Phase 3 (`dssatutils` acquisition)** — do for the new-site/new-crop use cases; keep the
   standalone fetch as fallback. Medium value, optional dep.
4. **Phase 4 (templates)** — do alongside the new-crop work; aligns scaffolding with the shared
   template dir.
5. **Output-mode B / engine daily parsers** — **defer**; only if the engine grows daily +
   sim-vs-meas parsing. Not worth coupling the calibrator's objective to a summary schema now.

**Net:** re-base the *plumbing* (execution → `dssatengine`, acquisition → `dssatutils`,
templates → `dssat_templates/`, inputs → `DSSAT48`) behind feature flags and optional extras;
**retain** the calibration-specific perturbation, daily/Evaluate parsing, objective and engines,
which no shared layer currently serves.

## 10. Implementation status (2026-06-18)

- `dssatengine` promoted public Python/R helpers: `run_dssat`,
  `write_dssbatch`, `write_dssbatch_sequence`, and treatment normalization
  (`normalize_treatment_list` in both languages). Private Python names remain
  aliases for older consumers. Version metadata is now `0.3.0`.
- `dssatcalibrator` added `execution.backend: native | dssatengine`; the shared
  backend delegates batch writing, treatment normalization, and DSSAT execution
  while keeping `dssat_io.py` PlantGro/Evaluate parsing local.
- `dssatcalibrator` formalized the `DSSAT48` install layout through
  `resolve_dssat_paths`.
- `weather.provider: dssatutils` and `soil.provider: dssatutils` provide the
  new-site acquisition path through optional `[acquire]` dependencies. The
  existing `file` and `nasa_power` paths remain fallbacks.
- `templates.template_dir` / `DSSAT_TEMPLATE_DIR` now resolve the shared
  `DSSAT_Gridded_Run_Tutorial/dssat_templates` convention, and crop scaffolding
  writes there by default when no `--out-dir` is given.
- Deferred by design: moving daily `PlantGro.OUT` and `Evaluate.OUT` parsers into
  `dssatengine`; the shared engine still exposes a summary-CSV contract.
