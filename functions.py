import numpy as np
import xarray as xr
import regionmask
import geopandas as gp
from importlib.resources import files
from functools import reduce
import pandas as pd

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

# ── CMIP6 Catalog and Model Discovery ─────────────────────────────────────────

def to_pystr_list(arr):
    """Convert a numpy array to a plain Python list of strings."""
    return arr.astype(str).tolist()

def preferred_load_list(cat_cloud, cat_esgf):
    """Split models into cloud-preferred and ESGF lists based on catalog availability."""
    cloud_models = reduce(np.union1d, [cat_cloud[t].df.source_id.unique() for t in ['tgt', 'awgt', 'aswgt']])
    esgf_models = reduce(np.union1d, [cat_esgf[t].df.source_id.unique() for t in ['tgt', 'awgt', 'aswgt']])
    all_models = np.union1d(cloud_models, esgf_models)
    cloud_subset = np.setdiff1d(all_models, esgf_models)
    return {
        'cloud': to_pystr_list(cloud_subset),
        'esgf': to_pystr_list(esgf_models),
    }

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

def _ann_Gt_old(da, member_dim="member_id", time_dim="time"):
    """Annual total in Gt yr⁻¹ (input: 10³ Gt month⁻¹). Climatology-first version, kept for comparison."""
    if da.sizes.get(member_dim, 0) == 0 or da.sizes.get(time_dim, 0) == 0:
        return 0.0
    clim = da.groupby(f"{time_dim}.month").mean(time_dim)
    return float(clim.sum("month").mean(member_dim).values) * 1e3

def _ann_Gt(da, member_dim="member_id", time_dim="time"):
    """Annual total in Gt yr⁻¹ (input: 10³ Gt month⁻¹). Sum-per-year, then mean over years."""
    if da.sizes.get(member_dim, 0) == 0 or da.sizes.get(time_dim, 0) == 0:
        return 0.0
    annual = da.groupby(f"{time_dim}.year").sum(time_dim)
    return float(annual.mean("year").mean(member_dim).values) * 1e3

def _t_safe(budget, key, s):
    """Annual Gt yr⁻¹ for a flux, or 0.0 if key is absent from budget."""
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
    inflows  = [(v, lb, col) for v, lb, col in inflows  if v > 0]
    outflows = [(v, lb, col) for v, lb, col in outflows if v > 0]
    imbalance = sum(v for v, *_ in inflows) - sum(v for v, *_ in outflows)
    total_throughput = sum(v for v, *_ in inflows) + sum(v for v, *_ in outflows)
    if total_throughput > 0 and abs(imbalance) / total_throughput > 0:
        dummy = (abs(imbalance), "Residual", "#aaaaaa")
        if imbalance > 0:
            outflows.append(dummy)
        else:
            inflows.append(dummy)
    return inflows, outflows

def _plotly_sankey_fig(res_label, res_color, inflows, outflows,
                       model_label, budget_type, title=None,
                       inflow_group=None, outflow_group=None, scale_ref=None,
                       show_values=True, show_residual=True, subtitle=None,
                       hide_terminal_nodes=False, show_percent=False):
    """
    Build a Plotly Sankey figure for a mass budget.

    inflows / outflows: list of (flow_Gt, label, color)

    inflow_group / outflow_group: optional dicts
      {"members": {label, ...}, "color": hex, "name": optional str}
      Routes the listed flows through an intermediate grouping node.

    scale_ref: optional (reservoir_total, n_busiest) tuple for shared figure
      height scaling across a set of Sankey figures.

    show_values: if False, omit Gt yr⁻¹ numbers from node labels.
    show_residual: if False, strip the Residual balancing node from the diagram.
    show_percent: if True, add each node's share of its side's total (inflows
      as % of total sources, outflows as % of total sinks) in parentheses.
    subtitle: optional string shown as a small annotation in the top-left corner.
    """
    import plotly.graph_objects as go

    inflows, outflows = _balance_flows(inflows, outflows)
    if not show_residual:
        inflows  = [(v, lb, col) for v, lb, col in inflows  if lb != "Residual"]
        outflows = [(v, lb, col) for v, lb, col in outflows if lb != "Residual"]
    n_in  = len(inflows)
    n_out = len(outflows)

    total_in  = sum(v for v, *_ in inflows)  or 1.0
    total_out = sum(v for v, *_ in outflows) or 1.0

    def _node_label(lb, v, total=None):
        parts = []
        if show_values:
            parts.append(f"{_fmt(v)} Gt yr⁻¹")
        if show_percent and total is not None:
            parts.append(f"{100 * v / total:.0f}%")
        if not parts:
            return lb
        return f"{lb} ({', '.join(parts)})" if lb else f"({', '.join(parts)})"

    node_labels = [_node_label(lb, v, total_in) for v, lb, _ in inflows]
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
            group_label = _node_label(group_name or "", combined_total) if show_values else (group_name or "")
            node_labels.append(group_label)
            node_colors.append(inflow_group["color"])
            node_x.append(0.25)
            node_y.append(combined_y)

    res_idx = len(node_labels)
    node_labels.append(res_label)
    node_colors.append(res_color)
    node_x.append(0.50)
    node_y.append(0.5)

    node_labels += [_node_label(lb, v, total_out) for v, lb, _ in outflows]
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
            out_group_label = _node_label(out_group_name or "", combined_out_total) if show_values else (out_group_name or "")
            node_labels.append(out_group_label)
            node_colors.append(outflow_group["color"])
            node_x.append(0.7)
            node_y.append(combined_out_y)

    if hide_terminal_nodes:
        transparent = "rgba(0,0,0,0)"
        for i in range(n_in):
            node_colors[i] = transparent
        for j in range(n_out):
            node_colors[res_idx + 1 + j] = transparent

    sources, targets, values, lcolors = [], [], [], []
    for i, (v, lb, col) in enumerate(inflows):
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
    for j, (v, lb, col) in enumerate(outflows):
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

    NODE_THICKNESS = 14
    NODE_PAD = 22
    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=NODE_PAD,
            thickness=NODE_THICKNESS,
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

    if title is not None:
        title_text = title
    else:
        model_id = "MME" if model_label == "Multi-model mean" else model_label
        prefix = f"{subtitle} " if subtitle else ""
        title_text = f"{prefix}{budget_type} Budget ({model_id})"

    fig.update_layout(
        title=dict(
            text=f"<b>{title_text}</b>",
            font=dict(size=14, color="#222222"),
        ),
        font=dict(size=11, color="#333333", family="Arial, sans-serif"),
        height=height,
        width=950,
        paper_bgcolor="white",
        margin=dict(l=20, r=20, t=60, b=120),
        annotations=[
            dict(x=0.47, y=1.06, xref="paper", yref="paper", text="sources",
                 showarrow=False, xanchor="right",
                 font=dict(size=11, color="#888888")),
            dict(x=0.53, y=1.06, xref="paper", yref="paper", text="sinks",
                 showarrow=False, xanchor="left",
                 font=dict(size=11, color="#888888")),
        ],
    )
    return fig


# ── Sankey: ice mass budget ────────────────────────────────────────────────────

def _ice_flows(ice_budget, model=None):
    """Build the (inflows, outflows) lists for an ice mass budget."""
    s = lambda da: _sel(da, model)
    bg_v = _t_safe(ice_budget, "basal growth",  s)
    fr_v = _t_safe(ice_budget, "frazil",        s)
    si_v = _t_safe(ice_budget, "snowice",       s)
    tm_v = _t_safe(ice_budget, "top melt",      s)
    bm_v = _t_safe(ice_budget, "basal melt",    s)
    lm_v = _t_safe(ice_budget, "lateral melt",  s)
    es_v = _t_safe(ice_budget, "evapsubl",      s)
    dy_v = _t_safe(ice_budget, "dynamics",      s)
    dyn_in  = max(0.0,  dy_v)
    dyn_out = max(0.0, -dy_v)

    inflows = [
        (abs(si_v), "Snow→ice",    _C["snow2ice_g"]),
        (abs(fr_v), "Open water ice production",      _C["frazil"]),
        (abs(bg_v), "Basal growth", _C["basal_g"]),
        (dyn_in,    "Dyn. import",  _C["dyn_in"]),
    ]
    outflows = [
        (abs(es_v), "Evap/subl",   _C["evapsubl"]),
        (abs(tm_v), "Top melt",     _C["top_melt"]),
        (dyn_out,   "Dyn. export",  _C["dyn_out"]),
        (abs(lm_v), "Lateral melt", _C["lat_melt"]),
        (abs(bm_v), "Basal melt",   _C["basal_melt"]),
    ]
    return inflows, outflows

def _ice_reservoir_total(ice_budget, model=None):
    """Return (reservoir_total, n_busiest) for an ice budget (used for shared scale refs)."""
    inflows, outflows = _balance_flows(*_ice_flows(ice_budget, model))
    return sum(v for v, *_ in inflows), max(len(inflows), len(outflows))

def make_ice_sankey_plotly(ice_budget, model=None, res_label=None, title=None,
                           inflow_group_name=None, outflow_group_name=None, scale_ref=None,
                           show_values=True, show_residual=True, subtitle=None, res_color=None,
                           hide_terminal_nodes=False, show_percent=False):
    """Ice mass budget as an interactive Plotly Sankey."""
    import plotly.graph_objects as go
    s = lambda da: _sel(da, model)
    if model is not None and s(ice_budget["basal growth"]).sizes["member_id"] == 0:
        return go.Figure()

    inflows, outflows = _ice_flows(ice_budget, model)
    inflow_group  = {"members": {"Open water ice production", "Basal growth"}, "color": _C["ice_growth"], "name": inflow_group_name}
    outflow_group = {"members": {"Basal melt", "Lateral melt"}, "color": _C["ice_melt"], "name": outflow_group_name}
    label = model if model is not None else "Multi-model mean"
    return _plotly_sankey_fig(res_label or "Sea Ice", res_color or _C["ice_res"],
                              inflows, outflows, label, "Sea Ice Mass", title=title,
                              inflow_group=inflow_group, outflow_group=outflow_group, scale_ref=scale_ref,
                              show_values=show_values, show_residual=show_residual, subtitle=subtitle,
                              hide_terminal_nodes=hide_terminal_nodes, show_percent=show_percent)


# ── Sankey: snow mass budget ───────────────────────────────────────────────────

def _snow_flows(snow_budget, model=None):
    """Build the (inflows, outflows) lists for a snow mass budget."""
    s = lambda da: _sel(da, model)
    sf_v = _t_safe(snow_budget, "snowfall",    s)
    sm_v = _t_safe(snow_budget, "snowmelt",    s)
    s2_v = _t_safe(snow_budget, "snow to ice", s)
    es_v = _t_safe(snow_budget, "evapsubl",    s)
    dy_v = _t_safe(snow_budget, "dynamics",    s)
    wd_v = _t_safe(snow_budget, "wind drift",  s)
    dyn_in   = max(0.0,  dy_v)
    dyn_out  = max(0.0, -dy_v)
    wind_in  = max(0.0,  wd_v)
    wind_out = max(0.0, -wd_v)

    inflows = [
        (abs(sf_v), "Snowfall",    _C["snowfall"]),
        (dyn_in,    "Dyn. import", _C["dyn_in"]),
        (wind_in,   "Wind import", _C["wind_in"]),
    ]
    outflows = [
        (abs(es_v),  "Evap/subl",  _C["evapsubl"]),
        (wind_out,   "Wind export", _C["wind_out"]),
        (abs(sm_v),  "Snowmelt",    _C["snowmelt"]),
        (abs(s2_v),  "Snow→ice",   _C["snow2ice_s"]),
        (dyn_out,    "Dyn. export", _C["dyn_out"]),
    ]
    return inflows, outflows

def _snow_reservoir_total(snow_budget, model=None):
    """Return (reservoir_total, n_busiest) for a snow budget (used for shared scale refs)."""
    inflows, outflows = _balance_flows(*_snow_flows(snow_budget, model))
    return sum(v for v, *_ in inflows), max(len(inflows), len(outflows))

def make_snow_sankey_plotly(snow_budget, model=None, res_label=None, title=None, scale_ref=None,
                            show_values=True, show_residual=True, subtitle=None, res_color=None,
                            hide_terminal_nodes=False, show_percent=False):
    """Snow mass budget as an interactive Plotly Sankey."""
    import plotly.graph_objects as go
    s = lambda da: _sel(da, model)
    if model is not None and s(snow_budget["snowfall"]).sizes["member_id"] == 0:
        return go.Figure()

    inflows, outflows = _snow_flows(snow_budget, model)
    label = model if model is not None else "Multi-model mean"
    return _plotly_sankey_fig(res_label or "Snow", res_color or _C["snow_res"],
                              inflows, outflows, label, "Snow on Sea Ice Mass", title=title,
                              scale_ref=scale_ref, show_values=show_values,
                              show_residual=show_residual, subtitle=subtitle,
                              hide_terminal_nodes=hide_terminal_nodes, show_percent=show_percent)

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
