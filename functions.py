import numpy as np
import xarray as xr
import regionmask
import geopandas as gp
from importlib.resources import files

def _data_path(file):
    return files('files').joinpath(file)

NH_seaice_regions = gp.read_file(_data_path('NSIDC-0780_SeaIceRegions_NH_v1.0.shp'))
SH_seaice_regions = gp.read_file(_data_path('NSIDC-0780_SeaIceRegions_SH-NASA_v1.0.shp'))

# ── Data loading helpers ───────────────────────────────────────────────────────

def preprocess_snap_to_month_start(ds):
    """Snap time coordinates to the first of each month for consistent merging."""
    new_time = ds.indexes['time'].to_period('M').to_timestamp()
    ds = ds.assign_coords(time=new_time)
    return ds


# ── Model selection utilities (from functions.py) ──────────────────────────────

def sel_model(ds, sid='CESM2'):
    """Select ensemble members for a single model by source_id prefix."""
    subset = ds.sel(member_id=ds.member_id.str.split('split', '_').sel(split=0) == sid)
    return subset

def drop_sel(ds: xr.Dataset, sids=('CESM2-LE',)):
    """Drop models whose member_id prefix matches any of the given sids."""
    if isinstance(sids, (str, bytes)):
        sids = [sids]
    prefix = ds['member_id'].astype(str).str.replace(r'_.+', '', regex=True)
    keep = ~prefix.isin(sids)
    if 'member_id' in ds.dims:
        return ds.isel(member_id=keep)
    else:
        members_to_keep = ds['member_id'].where(keep, drop=True)
        return ds.sel(member_id=members_to_keep)


# ── Sankey: seasonal ordering ──────────────────────────────────────────────────

START_MONTH = 3   # 3 = Mar (SH ice-growth onset);  7 = Jul (mid-winter)
_SH_MONTHS  = [(START_MONTH - 1 + i) % 12 + 1 for i in range(12)]
_MONTH_ABBR = {7:"Jul",8:"Aug",9:"Sep",10:"Oct",11:"Nov",12:"Dec",
               1:"Jan",2:"Feb",3:"Mar",4:"Apr",5:"May",6:"Jun"}
_SH_RANK    = {m: i for i, m in enumerate(_SH_MONTHS)}


# ── Sankey: model list and colour palette ──────────────────────────────────────

SANKEY_MODELS = ["ACCESS-CM2", "HadGEM3-GC31-LL", "UKESM1-0-LL",
                 "CESM2-LE",   "NorESM2-LM",       "NorESM2-MM"]

# Blue = thermodynamic gains   Red/orange = thermodynamic losses   Gray = dynamics
_C = {
    "ice_res":    "#0D47A1",
    "snow_res":   "#0277BD",
    "basal_g":    "#1565C0",
    "frazil":     "#1E88E5",
    "snow2ice_g": "#42A5F5",
    "snowfall":   "#0288D1",
    "top_melt":   "#B71C1C",
    "basal_melt": "#E53935",
    "lat_melt":   "#EF6C00",
    "evapsubl":   "#FF8F00",
    "snowmelt":   "#C62828",
    "snow2ice_s": "#FB8C00",
    "dyn_in":     "#616161",
    "dyn_out":    "#424242",
    "wind_in":    "#78909C",
    "wind_out":   "#546E7A",
    "ice_growth": "#1976D2",
    "ice_melt":   "#D84315",
}

_MIN_TRACE_Gt = 200   # Gt yr⁻¹ — below this: faint ribbon


# ── Sankey: flow computation helpers ──────────────────────────────────────────

def _sel(da, model):
    return sel_model(da, model) if model is not None else da

def _ann_Gt(da, member_dim="member_id", time_dim="time"):
    """Annual total in Gt yr⁻¹ (input: 10³ Gt month⁻¹)."""
    if da.sizes.get(member_dim, 0) == 0 or da.sizes.get(time_dim, 0) == 0:
        return 0.0
    clim = da.groupby(f"{time_dim}.month").mean(time_dim)
    return float(clim.sum("month").mean(member_dim).values) * 1e3

def _peak_sh(da, use_abs=False, member_dim="member_id", time_dim="time"):
    """Peak-month rank relative to START_MONTH and its abbreviation."""
    if da.sizes.get(member_dim, 0) == 0 or da.sizes.get(time_dim, 0) == 0:
        return 0, "N/A"
    clim = da.groupby(f"{time_dim}.month").mean(time_dim).mean(member_dim)
    arr  = np.abs(clim) if use_abs else clim
    pm   = int(arr.idxmax("month").values)
    return _SH_RANK[pm], _MONTH_ABBR[pm]

def _t(da, use_abs=False):
    v = _ann_Gt(da)
    r, pk = _peak_sh(da, use_abs=use_abs)
    return v, r, pk

def _t_safe(budget, key, s, use_abs=False):
    """Like _t but returns (0.0, 0, 'N/A') if key is absent from budget."""
    da = budget.get(key)
    if da is None:
        return 0.0, 0, "N/A"
    return _t(s(da), use_abs=use_abs)

def _ann_safe(budget, key, s):
    """Like _ann_Gt but returns 0.0 if key is absent from budget."""
    da = budget.get(key)
    return 0.0 if da is None else _ann_Gt(s(da))

def _fmt(v_Gt):
    v = abs(v_Gt)
    return f"{v/1e3:.1f}×10³" if v >= 500 else f"{v:.0f}"


# ── Sankey: Plotly rendering ───────────────────────────────────────────────────

def _hex_to_rgba(hex_color, alpha=0.50):
    h = hex_color.lstrip('#')
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"

def _balance_flows(inflows, outflows):
    """Filter zero flows and append a 'Residual' to whichever side is smaller."""
    inflows  = [(v, r, pk, lb, col) for v, r, pk, lb, col in inflows  if v > 0]
    outflows = [(v, r, pk, lb, col) for v, r, pk, lb, col in outflows if v > 0]
    imbalance = sum(v for v, *_ in inflows) - sum(v for v, *_ in outflows)
    total_throughput = sum(v for v, *_ in inflows) + sum(v for v, *_ in outflows)
    if total_throughput > 0 and abs(imbalance) / total_throughput > 0:
        dummy = (abs(imbalance), 99, "—", "Residual", "#aaaaaa")
        if imbalance > 0:
            outflows.append(dummy)
        else:
            inflows.append(dummy)
    return inflows, outflows

def _plotly_sankey_fig(res_label, res_color, inflows, outflows,
                       model_label, budget_type, title=None,
                       inflow_group=None, outflow_group=None, scale_ref=None):
    """
    Build a Plotly Sankey figure for a mass budget.

    inflows / outflows: list of (flow_Gt, peak_rank, peak_month, label, color)

    inflow_group / outflow_group: optional dicts
      {"members": {label, ...}, "color": hex, "name": optional str}
      Routes the listed flows through an intermediate grouping node.

    scale_ref: optional (reservoir_total, n_busiest) tuple for shared figure
      height scaling across a set of Sankey figures.
    """
    import plotly.graph_objects as go

    inflows, outflows = _balance_flows(inflows, outflows)
    n_in  = len(inflows)
    n_out = len(outflows)

    node_labels = [f"{lb} ({_fmt(v)} Gt yr⁻¹)" for v, _, _, lb, _ in inflows]
    node_colors = [col for *_, col in inflows]
    node_x = [0.01] * n_in

    def _flow_ys(terms):
        if not terms:
            return [0.5]
        flows = [v for v, *_ in terms]
        total = sum(flows) or 1.0
        pad = 0.02
        scale = 1.0 - 2 * pad
        cumulative = 0.0
        ys = []
        for f in flows:
            ys.append(pad + scale * (cumulative + f / 2) / total)
            cumulative += f
        return ys

    in_ys  = _flow_ys(inflows)
    out_ys = _flow_ys(outflows)
    node_y = list(in_ys)

    group_idx = None
    if inflow_group is not None:
        members = inflow_group["members"]
        grouped = [i for i, (*_, lb, _) in enumerate(inflows) if lb in members]
        if grouped:
            combined_total = sum(inflows[i][0] for i in grouped)
            combined_y = sum(node_y[i] * inflows[i][0] for i in grouped) / combined_total
            group_idx = len(node_labels)
            group_name = inflow_group.get("name")
            group_label = f"{group_name} ({_fmt(combined_total)} Gt yr⁻¹)" if group_name else f"{_fmt(combined_total)} Gt yr⁻¹"
            node_labels.append(group_label)
            node_colors.append(inflow_group["color"])
            node_x.append(0.25)
            node_y.append(combined_y)

    res_idx = len(node_labels)
    node_labels.append(res_label)
    node_colors.append(res_color)
    node_x.append(0.50)
    node_y.append(0.5)

    node_labels += [f"{lb} ({_fmt(v)} Gt yr⁻¹)" for v, _, _, lb, _ in outflows]
    node_colors += [col for *_, col in outflows]
    node_x += [0.99] * n_out
    node_y += out_ys

    out_group_idx = None
    if outflow_group is not None:
        out_members = outflow_group["members"]
        out_grouped = [j for j, (*_, lb, _) in enumerate(outflows) if lb in out_members]
        if out_grouped:
            combined_out_total = sum(outflows[j][0] for j in out_grouped)
            combined_out_y = sum(out_ys[j] * outflows[j][0] for j in out_grouped) / combined_out_total
            out_group_idx = len(node_labels)
            out_group_name = outflow_group.get("name")
            out_group_label = f"{out_group_name} ({_fmt(combined_out_total)} Gt yr⁻¹)" if out_group_name else f"{_fmt(combined_out_total)} Gt yr⁻¹"
            node_labels.append(out_group_label)
            node_colors.append(outflow_group["color"])
            node_x.append(0.7)
            node_y.append(combined_out_y)

    sources, targets, values, lcolors = [], [], [], []
    for i, (v, r, pk, lb, col) in enumerate(inflows):
        target = group_idx if (group_idx is not None and lb in inflow_group["members"]) else res_idx
        sources.append(i)
        targets.append(target)
        values.append(v)
        lcolors.append(_hex_to_rgba(col, 0.52))

    if group_idx is not None:
        sources.append(group_idx)
        targets.append(res_idx)
        values.append(combined_total)
        lcolors.append(_hex_to_rgba(inflow_group["color"], 0.52))

    out_start = res_idx + 1
    for j, (v, r, pk, lb, col) in enumerate(outflows):
        source = out_group_idx if (out_group_idx is not None and lb in outflow_group["members"]) else res_idx
        sources.append(source)
        targets.append(out_start + j)
        values.append(v)
        lcolors.append(_hex_to_rgba(col, 0.52))

    if out_group_idx is not None:
        sources.append(res_idx)
        targets.append(out_group_idx)
        values.append(combined_out_total)
        lcolors.append(_hex_to_rgba(outflow_group["color"], 0.52))

    hover = [f"{v:,.0f} Gt yr⁻¹" for v in values]

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=22,
            thickness=22,
            line=dict(color="white", width=0),
            label=node_labels,
            color=node_colors,
            x=node_x,
            y=node_y,
            hovertemplate="%{label}<extra></extra>",
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=lcolors,
            customdata=hover,
            hovertemplate="%{customdata}<extra></extra>",
        ),
    ))

    NODE_PAD = 22
    MARGIN_T, MARGIN_B = 60, 120
    PLOT_AREA_FULL = 480 - MARGIN_T - MARGIN_B  # 300

    n_busiest = max(n_in, n_out)
    if scale_ref is not None and scale_ref[0] > 0:
        ref_total, ref_n_busiest = scale_ref
        reservoir_total = sum(v for v, *_ in inflows)
        gaps_this = NODE_PAD * max(0, n_busiest - 1)
        gaps_ref  = NODE_PAD * max(0, ref_n_busiest - 1)
        available = (PLOT_AREA_FULL - gaps_ref) * (reservoir_total / ref_total)
        plot_area = max(50, gaps_this + available)
    else:
        plot_area = PLOT_AREA_FULL
    height = MARGIN_T + MARGIN_B + plot_area

    fig.update_layout(
        title=dict(
            text=f"<b>{title if title is not None else f'{budget_type} Budget: {model_label}'}</b>",
            font=dict(size=14, color="#222222"),
        ),
        font=dict(size=11, color="#333333", family="Arial, sans-serif"),
        height=height,
        width=950,
        paper_bgcolor="white",
        margin=dict(l=20, r=20, t=60, b=120),
    )
    return fig


# ── Sankey: ice mass budget ────────────────────────────────────────────────────

def _ice_flows(ice_budget, model=None):
    """Build the (inflows, outflows) lists for an ice mass budget."""
    s = lambda da: _sel(da, model)
    bg_v, bg_r, bg_pk = _t_safe(ice_budget, "basal growth",  s)
    fr_v, fr_r, fr_pk = _t_safe(ice_budget, "frazil",        s)
    si_v, si_r, si_pk = _t_safe(ice_budget, "snowice",       s)
    tm_v, tm_r, tm_pk = _t_safe(ice_budget, "top melt",      s, use_abs=True)
    bm_v, bm_r, bm_pk = _t_safe(ice_budget, "basal melt",   s, use_abs=True)
    lm_v, lm_r, lm_pk = _t_safe(ice_budget, "lateral melt", s, use_abs=True)
    es_v, es_r, es_pk = _t_safe(ice_budget, "evapsubl",      s, use_abs=True)
    dy_v, dy_r, dy_pk = _t_safe(ice_budget, "dynamics",      s, use_abs=True)
    dyn_in  = max(0.0,  _ann_safe(ice_budget, "dynamics", s))
    dyn_out = max(0.0, -_ann_safe(ice_budget, "dynamics", s))

    inflows = [
        (abs(si_v), si_r, si_pk, "Snow→Ice",    _C["snow2ice_g"]),
        (abs(fr_v), fr_r, fr_pk, "Frazil",      _C["frazil"]),
        (abs(bg_v), bg_r, bg_pk, "Basal Growth", _C["basal_g"]),
        (dyn_in,    dy_r, dy_pk, "Dyn. Import",  _C["dyn_in"]),
    ]
    outflows = [
        (abs(es_v), es_r, es_pk, "Evap/Subl",   _C["evapsubl"]),
        (abs(tm_v), tm_r, tm_pk, "Top Melt",     _C["top_melt"]),
        (dyn_out,   dy_r, dy_pk, "Dyn. Export",  _C["dyn_out"]),
        (abs(lm_v), lm_r, lm_pk, "Lateral Melt", _C["lat_melt"]),
        (abs(bm_v), bm_r, bm_pk, "Basal Melt",   _C["basal_melt"]),
    ]
    return inflows, outflows

def _ice_reservoir_total(ice_budget, model=None):
    """Return (reservoir_total, n_busiest) for an ice budget (used for shared scale refs)."""
    inflows, outflows = _balance_flows(*_ice_flows(ice_budget, model))
    return sum(v for v, *_ in inflows), max(len(inflows), len(outflows))

def make_ice_sankey_plotly(ice_budget, model=None, res_label=None, title=None,
                           inflow_group_name=None, outflow_group_name=None, scale_ref=None):
    """Ice mass budget as an interactive Plotly Sankey."""
    import plotly.graph_objects as go
    s = lambda da: _sel(da, model)
    if model is not None and s(ice_budget["basal growth"]).sizes["member_id"] == 0:
        return go.Figure()

    inflows, outflows = _ice_flows(ice_budget, model)
    inflow_group  = {"members": {"Frazil", "Basal Growth"}, "color": _C["ice_growth"], "name": inflow_group_name}
    outflow_group = {"members": {"Basal Melt", "Lateral Melt"}, "color": _C["ice_melt"], "name": outflow_group_name}
    label = model if model is not None else "Multi-model mean"
    return _plotly_sankey_fig(res_label or "Sea Ice", _C["ice_res"],
                              inflows, outflows, label, "Ice Mass", title=title,
                              inflow_group=inflow_group, outflow_group=outflow_group, scale_ref=scale_ref)


# ── Sankey: snow mass budget ───────────────────────────────────────────────────

def _snow_flows(snow_budget, model=None):
    """Build the (inflows, outflows) lists for a snow mass budget."""
    s = lambda da: _sel(da, model)
    sf_v, sf_r, sf_pk = _t_safe(snow_budget, "snowfall",    s)
    sm_v, sm_r, sm_pk = _t_safe(snow_budget, "snowmelt",    s, use_abs=True)
    s2_v, s2_r, s2_pk = _t_safe(snow_budget, "snow to ice", s, use_abs=True)
    es_v, es_r, es_pk = _t_safe(snow_budget, "evapsubl",    s, use_abs=True)
    dy_v, dy_r, dy_pk = _t_safe(snow_budget, "dynamics",    s, use_abs=True)
    wd_v, wd_r, wd_pk = _t_safe(snow_budget, "wind drift",  s, use_abs=True)
    dyn_in   = max(0.0,  _ann_safe(snow_budget, "dynamics",   s))
    dyn_out  = max(0.0, -_ann_safe(snow_budget, "dynamics",   s))
    wind_in  = max(0.0,  _ann_safe(snow_budget, "wind drift", s))
    wind_out = max(0.0, -_ann_safe(snow_budget, "wind drift", s))

    inflows = [
        (abs(sf_v), sf_r, sf_pk, "Snowfall",    _C["snowfall"]),
        (dyn_in,    dy_r, dy_pk, "Dyn. Import", _C["dyn_in"]),
        (wind_in,   wd_r, wd_pk, "Wind Import", _C["wind_in"]),
    ]
    outflows = [
        (abs(es_v),  es_r, es_pk, "Evap/Subl",  _C["evapsubl"]),
        (wind_out,   wd_r, wd_pk, "Wind Export", _C["wind_out"]),
        (abs(sm_v),  sm_r, sm_pk, "Snowmelt",    _C["snowmelt"]),
        (abs(s2_v),  s2_r, s2_pk, "Snow→Ice",   _C["snow2ice_s"]),
        (dyn_out,    dy_r, dy_pk, "Dyn. Export", _C["dyn_out"]),
    ]
    return inflows, outflows

def _snow_reservoir_total(snow_budget, model=None):
    """Return (reservoir_total, n_busiest) for a snow budget (used for shared scale refs)."""
    inflows, outflows = _balance_flows(*_snow_flows(snow_budget, model))
    return sum(v for v, *_ in inflows), max(len(inflows), len(outflows))

def make_snow_sankey_plotly(snow_budget, model=None, res_label=None, title=None, scale_ref=None):
    """Snow mass budget as an interactive Plotly Sankey."""
    import plotly.graph_objects as go
    s = lambda da: _sel(da, model)
    if model is not None and s(snow_budget["snowfall"]).sizes["member_id"] == 0:
        return go.Figure()

    inflows, outflows = _snow_flows(snow_budget, model)
    label = model if model is not None else "Multi-model mean"
    return _plotly_sankey_fig(res_label or "Snow", _C["snow_res"],
                              inflows, outflows, label, "Snow Mass", title=title,
                              scale_ref=scale_ref)

# ––––– Spatial Masking –––––––––––––––––

def region_mask(ds, region='Inner_Arctic'):
    """Mask a dataset to a named sea ice region."""
    ds = ds.copy()
    if 'lat' in ds.coords and 'y' in ds.dims:
        if region in ['Inner Arctic', 'IA']:
            region = 'Inner_Arctic'
        if region == 'Arctic':
            ds_subset = ds.where(ds.lat > 60)
        elif region == 'Antarctic':
            ds_subset = ds.where(ds.lat < -60)
        elif region == 'Inner_Arctic':
            df = NH_seaice_regions
            mask = regionmask.mask_geopandas(df, ds.lon, ds.lat, overlap=False)
            ds_subset = ds.where(mask.isin([0, 1, 2, 3, 4, 5, 6]))
        elif isinstance(region, dict):
            region_key = list(region.keys())[0]
            region_values = list(region.values())[0]
            if region_key == 'Arctic':
                df = NH_seaice_regions
            elif region_key == 'Antarctic':
                df = SH_seaice_regions
            else:
                raise ValueError(f"Unrecognized region dict key: {region_key}")
            mask = regionmask.mask_geopandas(df, ds.lon, ds.lat, overlap=False)
            ds_subset = ds.where(mask.isin(region_values))
        else:
            # If region doesn't match anything, just return unmasked
            ds_subset = ds
    else:
        print('No lat or y dim found: returning dataset')
        ds_subset = ds

    return ds_subset
