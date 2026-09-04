# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo does

Generates Sankey diagrams of sea ice and snow-on-ice mass budgets (sources: basal growth, frazil ice,
snowfall, snow→ice conversion; sinks: top/basal/lateral melt, evaporation/sublimation) from CMIP6 models,
the CESM2 Large Ensemble, and ERA5-driven reanalysis, area-integrated over the Arctic Ocean and Southern
Ocean for a 2015–2034 climatology (ERA5: 2000–2024). This is a research/analysis project built around
Jupyter notebooks and two shared Python modules — there is no application, build step, test suite, or CI.

## Environment setup

JupyterHub with existing scientific stack — just add the missing packages:
```bash
pip install cf-xarray xesmf esgf-pyclient intake-esgf globus-sdk cmocean kaleido
```

Fresh conda/mamba env (preferred for `xesmf`/`cartopy` binary deps):
```bash
mamba env create -f environment.yml
conda activate sea-ice-mass-balance
python -m ipykernel install --user --name sea-ice-mass-balance --display-name "sea-ice-mass-balance"
```

Fresh pip/uv env: `uv venv --python 3.12 && uv pip install -r requirements.txt` (same deps as `environment.yml`,
kept in sync manually — update both when changing dependencies).

There are no lint/test/build commands in this repo. "Running" the project means executing notebook cells.

## Architecture: two-stage pipeline

**Stage 1 — `model_load.ipynb`** (uses `load.py`): Loads gridded CMIP6 budget variables from ESGF/Pangeo
cloud catalogs (via `intake-esgf`, falling back to direct ESGF downloads with retry logic) and CESM2-LE
from local storage or THREDDS OPeNDAP. Applies all model-specific sign/unit/area-basis corrections (see
README "Model-Specific Corrections" table — these are load-time patches for known bugs/quirks in specific
models' published output, not general-purpose logic). Area-integrates fluxes over the Southern Ocean,
Weddell Sea, and Arctic Ocean using `regionmask`/NSIDC-0780 region shapefiles. Saves one `.nc` file per
variable per model to `save_path`. Intended to be run one CMIP6 model at a time (set `models`).

**Stage 2 — `cmip6_sankey.ipynb`** (uses `functions.py`): Reads the `.nc` files written by stage 1
from `melt_path` (must match stage 1's `save_path`), merges CMIP6 and CESM2-LE data into per-region budget
dicts (`ice_budget_*`, `snow_budget_*`, now including an Inner Arctic region: `ice_budget_IA`), computes a
shared `scale_ref` height scale so figures are visually comparable, and calls `make_ice_sankey_plotly` /
`make_snow_sankey_plotly` once per model in `SANKEY_MODELS` plus once with `model=None` for the multi-model
ensemble (MME) mean. `SANKEY_MODELS` is defined in the notebook itself (not `functions.py`) and is
conditional on an `add_CESM2_LE` toggle — CESM2-LE is off by default since CESM2 + CESM2-WACCM already
cover that model family and including all three would over-weight it in the MME. Figures are written to
`figures/` as `.png`/`.pdf` (optionally `.svg`), named `{ice,snow}_sankey_{AO,SO,IA}_{model}_{year-range}.{ext}`.

**Independent — `era5_sankey.ipynb`**: builds a comparison Sankey from an ERA5-forced NEMO-SI3
ocean–sea-ice simulation (eORCA025, 2000–2024, data from Benjamin Richaud), used to sanity-check the CMIP6
Sankeys against a forced-reanalysis baseline. Does not go through the `model_load.ipynb`/`load.py` pipeline.

Both notebooks import from the shared modules rather than duplicating logic — when editing budget-term
definitions, unit conversions, or corrections, change `load.py`/`functions.py`, not notebook cells.

### `load.py` — data acquisition and preprocessing
- `CMIP6` and `CESM` classes are the main entry points for pulling/opening a model's budget variables.
- Catalog/download helpers (`load_from_catalog`, `_download_and_open`, `_download_file_with_retries`,
  THREDDS catalog helpers) handle ESGF/cloud/THREDDS access with local caching and retry logic.
- Preprocessing helpers (`complete_preprocessing`, `sanitize_time`, `convert_time`, `add_corners`,
  `calc_areacello`) normalize grids/time coordinates across heterogeneous CMIP6 model output before it can
  be merged or regridded (uses `xmip` combined preprocessing under the hood).
- `NH_seaice_regions` / `SH_seaice_regions` (from `files/NSIDC-0780_SeaIceRegions_*.shp`) define region
  polygons used for masking/integration.

### `functions.py` — budget assembly and Sankey rendering
- `_C` is the shared color palette (blue = thermodynamic gains, red/orange = thermodynamic losses, gray =
  dynamics) — keep new flow types consistent with this convention rather than picking arbitrary colors.
  (The model roster, `SANKEY_MODELS`, now lives in `cmip6_sankey.ipynb`, not here.)
- `_ann_Gt` converts monthly flux data (input units: 10³ Gt/month) to an annual Gt/yr climatology by
  summing per year then averaging across years (preferred over the deprecated climatology-first
  `_ann_Gt_old`, kept only for comparison).
- `_ice_flows`/`_snow_flows` build the (value, label, color) flow lists per budget from a `budget` dict of
  DataArrays (keyed by term name, e.g. `'basal growth'`, `'snowfall'`); `_balance_flows` adds a "Residual"
  flow to whichever side (sources/sinks) is smaller so the diagram balances.
- `_plotly_sankey_fig` is the low-level Plotly Sankey builder; `make_ice_sankey_plotly`/
  `make_snow_sankey_plotly` wrap it for the two budget types. `scale_ref` is threaded through so a set of
  figures shares one height scale — always pass the notebook's precomputed `ICE_SANKEY_SCALE_REF`/
  `SNOW_SANKEY_SCALE_REF` rather than letting each figure size independently.
- `sel_model`/`drop_sel` filter a merged multi-model dataset by `member_id` prefix (e.g. to exclude a model
  from a figure, or to include/exclude CESM2-LE via `add_CESM2_LE`).
- `region_mask` applies the NSIDC-0780 region polygons for masking a dataset to a named region (default
  `'Inner_Arctic'`).
- `plot_ice_terms_climatology` (in `cmip6_sankey.ipynb`) is a diagnostic multipanel plot — one
  panel per ice budget term, one line per model's monthly climatology — used to spot-check corrections
  across models rather than to produce a paper figure.

## Models

Nine models are used by default (up from the original six): ACCESS-CM2, HadGEM3-GC31-LL, UKESM1-0-LL,
MRI-ESM2-0, CESM2-WACCM, NorESM2-LM, CNRM-CM6-1, CNRM-ESM2-1, IPSL-CM6A-LR. CESM2, NorESM2-MM, and
CNRM-CM6-1-HR have also been ingested but are excluded from the default set as resolution/model-family
duplicates of CESM2-WACCM, NorESM2-LM, and CNRM-CM6-1 respectively. CESM2-LE is supported but off by
default (`add_CESM2_LE=True` in `cmip6_sankey.ipynb` to include it) — with CESM2-WACCM already in
the default set, adding CESM2-LE too would over-weight that model family in the MME.

## Key data-modeling gotchas

- **Sign convention**: budget terms use losses-negative / gains-positive throughout. Several models report
  some loss terms as positive and require a `-1` correction applied in `model_load.ipynb` — see README
  "Model-Specific Corrections" before adding a new model or trusting raw values.
- **Area basis**: some models (ACCESS-CM2, and all its snow variables) report fluxes per sea-ice area
  rather than per grid-cell area and must be multiplied by `siconc/100` before area-integration.
- **Units**: `_ann_Gt` assumes input is 10³ Gt/month; NorESM2 `sndmasssnf` needs a divide-by-330
  (kg m⁻³ snow density) correction due to a known upstream CMOR bug
  ([NorESMhub/noresm2cmor#282](https://github.com/NorESMhub/noresm2cmor/issues/282)).
- **CNRM basal growth/top melt are approximated**: those models don't report the two terms separately, so
  `sidmassgrowthbot` and `sidmassmelttop` are reconstructed by summing the raw terms and splitting by sign
  — a cell with both growth and melt in the same month gets misattributed entirely to whichever dominates.
- **CESM2-WACCM has per-member archiving bugs**, not model-wide ones: member `r1i1p1f1` has several ice
  terms archived 1800× too small (rescaled ×1800) while every other member has `sidmassevapsubl` 1e6× too
  large (rescaled ÷1e6); on the snow side `r1i1p1f1`'s snowmelt is 1800× too small and `r3i1p1f1`'s
  snowfall is ~330× too large. Check which member you're looking at before trusting a raw value.
- **Some corrections are diagnosed, not confirmed by the modeling centers** — the ACCESS-CM2 snowmelt fix
  (density-unit mismatch, corrected ÷(1000/330)), the CNRM growth/melt split, and CESM2-WACCM's per-member
  rescaling are all inferred from output behavior. Treat them as best-effort when interpreting results.
- **`melt_path`/`save_path` must match** between `model_load.ipynb` and `cmip6_sankey.ipynb` — the
  second notebook only sees variables the first one actually wrote out.
- `data/` (raw/intermediate `.nc` outputs) is gitignored; only final figures and static reference files
  under `files/` are tracked.
