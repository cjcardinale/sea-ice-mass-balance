# Sea Ice Mass Balance Sankey Diagrams

**Lead:** Alek Petty | **Data visuals:** Chris Cardinale | **Co-I:** Maddie Smith

With thanks to the SIMIP community for helpful discussion.

## Environment Setup

### JupyterHub / existing conda environment

If you are working in a shared JupyterHub where the base environment already has most scientific packages (numpy, xarray, dask, geopandas, etc.), install the missing packages with pip:

```bash
pip install cf-xarray xesmf esgf-pyclient intake-esgf globus-sdk cmocean kaleido
```

### Fresh install: conda/mamba (recommended)

`xesmf` and `cartopy` have binary dependencies that conda resolves more reliably than pip. [Mamba](https://mamba.readthedocs.io) is a faster drop-in replacement for conda.

```bash
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate sea-ice-mass-balance
python -m ipykernel install --user --name sea-ice-mass-balance --display-name "sea-ice-mass-balance"
```

Then restart JupyterHub/JupyterLab and select the `sea-ice-mass-balance` kernel.

### Fresh install: pip/uv

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv venv --python 3.12
source .venv/bin/activate
uv pip install -r requirements.txt
python -m ipykernel install --user --name sea-ice-mass-balance --display-name "sea-ice-mass-balance"
```

## Overview

Sankey diagrams of sea ice and snow-on-ice mass budgets derived from CMIP6 models and the CESM2 Large Ensemble (CESM2-LE). Budget terms are area-integrated over the Southern Ocean and Arctic Ocean for a 2015–2034 climatological period, then displayed as flow diagrams where ribbon width is proportional to annual mass flux. Flux widths are normalized across both hemispheres so Arctic and Southern Ocean figures can be compared directly.

The diagrams show the primary sources (basal growth, frazil ice, snowfall, snow→ice conversion) and sinks (top melt, basal melt, lateral melt, evaporation/sublimation) of each reservoir, with thermodynamic growth terms grouped into an intermediate **Ocean Growth** node and thermodynamic melt terms into an **Ocean Melt** node.

### Ice Mass Budget
**Southern Ocean**
![SO Ice MME](figures/ice_sankey_SO_MME_2015-2034.png)

**Arctic Ocean**
![AO Ice MME](figures/ice_sankey_AO_MME_2015-2034.png)

### Snow Mass Budget
**Southern Ocean**
![SO Snow MME](figures/snow_sankey_SO_MME_2015-2034.png)

**Arctic Ocean**
![AO Snow MME](figures/snow_sankey_AO_MME_2015-2034.png)

## Models and Period

| Model | Institution | Notes |
|---|---|---|
| ACCESS-CM2 | CSIRO/ARCCSS | Requires area-basis, sign, and snowmelt corrections (see [Model-Specific Corrections](#model-specific-corrections)) |
| HadGEM3-GC31-LL | Met Office | Snow dynamics available; no snow wind drift or snow sublimation output |
| UKESM1-0-LL | Met Office / MOHC | Snow dynamics available; no snow wind drift or snow sublimation output |
| NorESM2-LM | NCC | Snow wind drift available; no snow dynamics or snow sublimation; requires sign/unit corrections |
| NorESM2-MM | NCC | Snow wind drift available; no snow dynamics or snow sublimation; requires sign/unit corrections |
| CESM2 | NCAR | Requires melt-term sign corrections |
| CESM2-WACCM | NCAR | Requires the same corrections as CESM2, plus per-member archiving fixes |
| CESM2-LE | NCAR | Large ensemble (100 members); **optional, off by default** — superseded by CESM2 + CESM2-WACCM in the default multi-model mean (set `add_CESM2_LE=True` in `sankey_figures_clean.ipynb` to include it) |
| MRI-ESM2-0 | MRI | No corrections needed; no snow wind drift output |
| CNRM-CM6-1 | CNRM-CERFACS | Basal growth/top melt are approximated (see caveats); no snow wind drift output |
| CNRM-CM6-1-HR | CNRM-CERFACS | Basal growth/top melt are approximated (see caveats); no snow wind drift output |
| CNRM-ESM2-1 | CNRM-CERFACS | Basal growth/top melt are approximated (see caveats); no snow wind drift output |

**Period:** 2015–2034 (20-year climatological mean across all available ensemble members)

## Workflow

Two notebooks run in sequence:

1. **`model_load.ipynb`** — Loads gridded CMIP6 budget variables from ESGF/Pangeo and CESM2-LE from local storage or THREDDS OPeNDAP. Applies all model-specific corrections (see below). Area-integrates fluxes over the Southern Ocean, Weddell Sea, and Arctic Ocean. Saves one `.nc` file per variable per model to `save_path`. Process one CMIP6 model at a time by setting `models`.

2. **`sankey_figures_clean.ipynb`** — Reads the `.nc` files from `melt_path` (must match `save_path`), merges CMIP6 and CESM2-LE data, and generates Sankey diagrams for each model and the multi-model ensemble mean. Figures are written to `figures/` as both `.png` and `.pdf` and optionally `.svg`.

Separately, **`era5_sankey.ipynb`** builds a comparison Sankey from an ERA5-forced NEMO-SI3 ocean–sea ice simulation (eORCA025, 2000–2024, data provided by Benjamin Richaud) — independent of the CMIP6 loading pipeline above, used to sanity-check the model Sankeys against a forced-reanalysis baseline.

## Model-Specific Corrections

Several models required corrections to standardize outputs to a common sign convention (losses negative, gains positive) and a common area basis (per grid-cell area) before area-integration. These are applied in `model_load.ipynb`.

### Ice Mass Budget

| Model | Correction |
|---|---|
| **ACCESS-CM2** | `sidmassth` and `sidmassdyn` are per grid-cell area (no scaling needed). All other ice flux variables are incorrectly reported per sea-ice area — multiplied by `siconc/100` to convert to per grid-cell area before integration. `sidmassevapsubl` is reported positive (mass loss) — multiplied by −1 to enforce the loss-negative sign convention. |
| **HadGEM3-GC31-LL, UKESM1-0-LL** | `sidmassevapsubl` is reported positive (mass loss) — multiplied by −1 to enforce the loss-negative sign convention. |
| **NorESM2-LM, NorESM2-MM, CESM2, CESM2-LE** | `sidmassmelttop`, `sidmassmeltbot`, and `sidmasslat` are reported positive (mass loss) — multiplied by −1 to enforce the loss-negative sign convention. |
| **CESM2-WACCM** | Same melt-term sign correction as CESM2. Per-member archiving issues, not model-wide: only member `r1i1p1f1` has `sidmassgrowthbot`, `sidmassgrowthwat`, `sidmasssi`, `sidmassmelttop`, `sidmassmeltbot`, and `sidmasslat` archived 1800× too small (consistent with CESM2's 30-minute ice–ocean coupling interval, though the root cause hasn't been confirmed) — that member is rescaled ×1800, others are left as archived. Conversely, `sidmassevapsubl` is 1e6× too large in every member *except* `r1i1p1f1` — those members are rescaled ÷1e6. `sidmassth` and `sidmassdyn` needed no correction in any member. |
| **MRI-ESM2-0** | No corrections needed — archived terms already follow the CF sign convention. |
| **CNRM-CM6-1, CNRM-CM6-1-HR, CNRM-ESM2-1** | `sidmassgrowthbot` (basal growth) and `sidmassmelttop` (top melt) are reconstructed from their sum: the two raw terms are added together, then the combined field is split by sign — positive values assigned to basal growth, negative values to top melt. This is an approximation: a grid cell with both growth and melt in the same month is misattributed entirely to whichever process dominates. All other ice terms need no correction. |

### Snow Mass Budget

| Model | Correction |
|---|---|
| **ACCESS-CM2** | All snow variables reported per sea-ice area — multiplied by `siconc/100`. Snowmelt (`sndmassmelt`) is anomalously large: CICE's cm/day → kg/s conversion should use snow density (330 kg m⁻³), but ACCESS-CM2 appears to have used freshwater density (1000 kg m⁻³) instead — corrected by dividing by `1000/330`. This diagnosis isn't confirmed in ACCESS-CM2's own documentation, but the corrected magnitudes now align with other models. `sndmasswindrif` and `sndmasssubl` are unavailable. |
| **HadGEM3-GC31-LL, UKESM1-0-LL** | `sndmasswindrif` and `sndmasssubl` are unavailable. |
| **NorESM2-LM, NorESM2-MM** | `sndmassmelt` is reported positive (mass loss) — multiplied by −1 to enforce the loss-negative sign convention. `sndmasssnf` values are 330× too large due to a unit conversion error in the published CMIP6 output — divided by snow density (330 kg m⁻³) to correct ([NorESMhub/noresm2cmor #282](https://github.com/NorESMhub/noresm2cmor/issues/282)). `sndmasssubl` and `sndmassdyn` are unavailable. |
| **CESM2** | `sndmassmelt` is reported positive (mass loss) — multiplied by −1, same as CESM2-LE. |
| **CESM2-WACCM** | Same snowmelt sign correction as CESM2. Per-member archiving issues, not model-wide: only member `r1i1p1f1` has archived snowmelt 1800× too small (same coupling-interval issue as the ice budget) — rescaled ×1800, other members left as archived. Separately, only member `r3i1p1f1` has archived snowfall ~330× too large (same snow-density-unit issue as the NorESM2 fix) — rescaled ÷330, other members' snowfall is archived correctly. `sndmasssubl` is available (unlike most models) but left unscaled since its magnitude hasn't been checked; snow-to-ice and wind drift are also not yet checked. |
| **CESM2-LE** | `sndmassmelt` is reported positive (mass loss) — multiplied by −1 to enforce the loss-negative sign convention. Snow→ice transfer (`sidmasssi`) is on the ice-budget side — converted to snow-budget units by multiplying by `-(330/917)` (snow/ice density ratio). |
| **MRI-ESM2-0, CNRM-CM6-1, CNRM-CM6-1-HR, CNRM-ESM2-1** | No corrections needed — archived terms already follow the CF sign convention. `sndmasswindrif` is unavailable. |

### Variable Availability by Model

Availability across the full CMIP6 archive (native grid, `SImon`, `ssp245`) as surfaced by the catalog search in `model_load.ipynb` (`availability_df`), plus CESM2-LE (loaded separately from NCAR storage, not through the ESGF/cloud catalog). ⭐ marks the models used in the Sankey figures (alphabetical order); CESM2-LE is supported but off by default (see [Known Issues and Caveats](#known-issues-and-caveats)).

**Data availability summary — 36 total models**

| Ice Variable | # Models with data |
|---|:---:|
| Top Melt | 30 |
| Basal Melt | 30 |
| Lateral Melt | 22 |
| Basal Growth | 28 |
| Frazil | 29 |
| Snow→Ice | 30 |
| Evap/Subl | 28 |
| Dynamics | 22 |

| Snow Variable | # Models with data |
|---|:---:|
| Melt | 25 |
| Snowfall | 33 |
| Snow→Ice | 15 |
| Evap/Subl | 10 |
| Wind Drift | 2 |
| Dynamics | 11 |

**Ice Mass Budget**

| Model | Top Melt | Basal Melt | Lateral Melt | Basal Growth | Frazil | Snow→Ice | Evap/Subl | Dynamics |
|---|---|---|---|---|---|---|---|---|
| ⭐ ACCESS-CM2 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ACCESS-ESM1-5 | — | — | — | — | — | — | ✓ | — |
| AWI-CM-1-1-MR | — | — | — | — | — | ✓ | ✓ | — |
| BCC-CSM2-MR | ✓ | ✓ | — | — | ✓ | ✓ | ✓ | — |
| CAMS-CSM1-0 | ✓ | ✓ | — | — | — | ✓ | ✓ | — |
| CAS-ESM2-0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| ⭐ CESM2 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ⭐ CESM2-LE | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ⭐ CESM2-WACCM | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| CIESM | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| CMCC-CM2-SR5 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| CMCC-ESM2 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ⭐ CNRM-CM6-1 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ⭐ CNRM-CM6-1-HR | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ⭐ CNRM-ESM2-1 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| EC-Earth3 | ✓ | ✓ | — | ✓ | ✓ | ✓ | ✓ | ✓ |
| EC-Earth3-CC | ✓ | ✓ | — | ✓ | ✓ | ✓ | ✓ | — |
| EC-Earth3-Veg | ✓ | ✓ | — | ✓ | ✓ | ✓ | ✓ | — |
| EC-Earth3-Veg-LR | ✓ | ✓ | — | ✓ | ✓ | ✓ | ✓ | — |
| FGOALS-f3-L | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| FGOALS-g3 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| FIO-ESM-2-0 | — | — | — | — | — | — | — | — |
| GISS-E2-1-G | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| GISS-E2-1-G-CC | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| GISS-E2-1-H | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| GISS-E2-2-G | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| ⭐ HadGEM3-GC31-LL | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| IPSL-CM6A-LR | ✓ | ✓ | — | ✓ | ✓ | ✓ | ✓ | ✓ |
| MPI-ESM1-2-HR | — | — | — | — | — | — | — | ✓ |
| MPI-ESM1-2-LR | — | — | — | — | — | — | — | ✓ |
| ⭐ MRI-ESM2-0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| NESM3 | — | — | — | — | — | — | — | — |
| ⭐ NorESM2-LM | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ⭐ NorESM2-MM | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| TaiESM1 | ✓ | ✓ | — | ✓ | ✓ | — | — | — |
| ⭐ UKESM1-0-LL | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

**Snow Mass Budget**

| Model | Melt | Snowfall | Snow→Ice | Evap/Subl | Wind Drift | Dynamics |
|---|---|---|---|---|---|---|
| ⭐ ACCESS-CM2 | ✓ | ✓ | ✓ | — | — | ✓ |
| ACCESS-ESM1-5 | — | ✓ | — | — | — | — |
| AWI-CM-1-1-MR | — | — | — | — | — | — |
| BCC-CSM2-MR | ✓ | ✓ | — | — | — | — |
| CAMS-CSM1-0 | — | ✓ | — | — | — | — |
| CAS-ESM2-0 | ✓ | ✓ | — | — | — | — |
| ⭐ CESM2 | ✓ | ✓ | — | — | — | — |
| ⭐ CESM2-LE | ✓ | ✓ | ✓ | ✓ | — | — |
| ⭐ CESM2-WACCM | ✓ | ✓ | — | ✓ | — | — |
| CIESM | ✓ | ✓ | — | — | — | — |
| CMCC-CM2-SR5 | ✓ | ✓ | — | ✓ | — | — |
| CMCC-ESM2 | ✓ | ✓ | — | ✓ | — | — |
| ⭐ CNRM-CM6-1 | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| ⭐ CNRM-CM6-1-HR | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| ⭐ CNRM-ESM2-1 | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| EC-Earth3 | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| EC-Earth3-CC | ✓ | ✓ | — | — | — | — |
| EC-Earth3-Veg | ✓ | ✓ | — | — | — | — |
| EC-Earth3-Veg-LR | ✓ | ✓ | — | — | — | — |
| FGOALS-f3-L | ✓ | ✓ | ✓ | — | — | — |
| FGOALS-g3 | ✓ | ✓ | ✓ | — | — | — |
| FIO-ESM-2-0 | — | ✓ | — | — | — | — |
| GISS-E2-1-G | — | ✓ | — | — | — | — |
| GISS-E2-1-G-CC | — | ✓ | — | — | — | — |
| GISS-E2-1-H | — | ✓ | — | — | — | — |
| GISS-E2-2-G | — | ✓ | — | — | — | — |
| ⭐ HadGEM3-GC31-LL | ✓ | ✓ | ✓ | — | — | ✓ |
| IPSL-CM6A-LR | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| MPI-ESM1-2-HR | — | ✓ | — | — | — | ✓ |
| MPI-ESM1-2-LR | — | ✓ | — | — | — | ✓ |
| ⭐ MRI-ESM2-0 | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| NESM3 | — | — | ✓ | — | — | — |
| ⭐ NorESM2-LM | ✓ | ✓ | ✓ | — | ✓ | — |
| ⭐ NorESM2-MM | ✓ | ✓ | ✓ | — | ✓ | — |
| TaiESM1 | ✓ | — | — | — | — | — |
| ⭐ UKESM1-0-LL | ✓ | ✓ | ✓ | — | — | ✓ |

**Note:** Snow-side snow→ice conversion (`sndmasssi`) is missing for many models, but can be derived from the ice-side snow-to-ice tendency (`sidmasssi`) via the snow/ice density ratio when the ice-side variable is available. This is how CESM2-LE's `sndmasssi` is computed (see [Model-Specific Corrections](#snow-mass-budget)).

## Known Issues and Caveats

- **Model results only.** These are not observationally constrained budgets.
- **Eleven models are used by default**, up from the original six. CESM2-LE is supported but off by default, since CESM2 and CESM2-WACCM are already in the default set and including all three would over-weight the CESM2 family in the multi-model mean.
- **Snow budgets are uncertain.** Inconsistent and incomplete snow flux outputs across models make the snow Sankeys less reliable than the ice Sankeys.
- **Some corrections are diagnosed rather than confirmed** — the ACCESS-CM2 snowmelt fix, the CNRM basal-growth/top-melt split, and CESM2-WACCM's per-member rescaling are all inferred from the archived output rather than documented by the modeling centers. See [Model-Specific Corrections](#model-specific-corrections) for what each one assumes and how confident it is.
- **Non-zero wind drift and dynamics at hemispheric scale.** Summing over a full hemisphere should produce near-zero dynamics and wind drift, but small non-zero values remain. This may reflect physical processes (e.g., snow lost to the ocean via wind or ice deformation at the ice edge) or model artifacts.
- **Evaporation/sublimation** on the snow side is available for CESM2-LE, CESM2-WACCM, MRI-ESM2-0, and the CNRM models; other models do not report it.
- **Wind drift** is only available from NorESM2 models.
- **Ice dynamics** on the snow side is only available from HadGEM3-GC31-LL and UKESM1-0-LL.
- **Frazil ice** calculations are inconsistent across models and require further investigation.

## Updates

**August 11, 2026**
- Expanded the default model set from six to eleven: added CESM2, CESM2-WACCM, MRI-ESM2-0, CNRM-CM6-1, CNRM-CM6-1-HR, and CNRM-ESM2-1. CESM2-LE remains supported but off by default.
- Documented each new model's corrections (or lack thereof) in [Model-Specific Corrections](#model-specific-corrections), including CESM2-WACCM's per-member archiving fixes and the CNRM basal-growth/top-melt approximation.
- Re-diagnosed the ACCESS-CM2 snowmelt correction as a snow-density-vs-freshwater-density unit mixup (divide by `1000/330`), replacing the earlier `/3.3` heuristic.
- Added `era5_sankey.ipynb`, an ERA5-forced NEMO-SI3 comparison Sankey (2000–2024).

**June 12, 2026**
- Flux bar widths normalized across hemispheres to enable direct Arctic/Southern Ocean comparisons.
- Added intermediate **Ocean Growth** and **Ocean Melt** grouping nodes to connect thermodynamically related terms and reflect differing model definitions in the literature.

## Extra plot ideas

- Follow the Keen et al. (2021) sea ice mass budget analysis and include an Inner Arctic Ocean domain approach that would allow for a dynamics/ice export contribution.
