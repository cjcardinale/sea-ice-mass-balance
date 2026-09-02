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

# ── Model selection utilities ───────────────────────────────────────────────────

def sel_model(ds, sid='CESM2'):
    """Select ensemble members for a single model by source_id prefix."""
    subset = ds.sel(member_id=ds.member_id.str.split('split', '_').sel(split=0) == sid)
    return subset

def drop_sel(ds: xr.Dataset, sids=('CESM2-LE',)):
    """Drop models whose member_id prefix matches any of the given sids."""
    if ds is None:
        return None
    if isinstance(sids, (str, bytes)):
        sids = [sids]
    prefix = ds['member_id'].astype(str).str.replace(r'_.+', '', regex=True)
    keep = ~prefix.isin(sids)
    if 'member_id' in ds.dims:
        return ds.isel(member_id=keep)
    else:
        members_to_keep = ds['member_id'].where(keep, drop=True)
        return ds.sel(member_id=members_to_keep)


# ── Sankey: color palette ──────────────────────────────────────

# Blue = thermodynamic gains   Red/orange = thermodynamic losses   Gray = dynamics
_C = {
    # Ocean-contact terms grouped into "Ocean Growth"/"Ocean Melt" in the ice Sankey
    # keep tight blue/red families anchored on the group node's own color.
    "ice_res":    "#0D47A1",
    "snow_res":   "#0277BD",
    "ice_growth": "#1976D2",   # "Ocean Growth" group node
    "basal_g":    "#1565C0",   # in group — blue family
    "frazil":     "#1E88E5",   # in group — blue family
    "ice_melt":   "#D84315",   # "Ocean Melt" group node
    "basal_melt": "#B71C1C",   # in group — red family
    "lat_melt":   "#E53935",   # in group — red family
    # Non-ocean-contact terms get their own families: snow→ice is an internal
    # ice/snow conversion (teal); top melt and evap/sublimation are atmosphere-driven
    # (gold/orange).
    "snow2ice_g": "#00897B",
    "top_melt":   "#E65100",
    "evapsubl":   "#EF6C00",
    # Snow-side terms (not part of any ice grouping node).
    "snowfall":   "#0288D1",
    "snowmelt":   "#C62828",
    "snow2ice_s": "#FB8C00",
    "wind_in":    "#78909C",
    "wind_out":   "#546E7A",
    # Dynamics (transport, not a source/sink) — shared by ice and snow.
    "dyn_in":     "#616161",
    "dyn_out":    "#424242",
}

_MIN_TRACE_Gt = 200   # Gt yr⁻¹ — below this: faint ribbon


# ── Sankey: flow computation helpers ──────────────────────────────────────────

def _sel(da, model):
    """Restrict to one model's members, or pass through unchanged if model is None."""
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

def _drop_placeholder_zero_members(da, member_dim="member_id", time_dim="time"):
    """Drop members whose data is identically zero across the whole time series —
    a real flux essentially never sums to exactly zero, so this is a placeholder
    (a model that doesn't archive the variable) rather than a genuine zero. Left
    in, it would silently dilute the multi-model mean for that term."""
    if da.sizes.get(member_dim, 0) == 0 or da.sizes.get(time_dim, 0) == 0:
        return da
    nonzero = (da != 0).any(time_dim).compute()
    return da.where(nonzero, drop=True)

def _t_safe(budget, key, s):
    """Annual Gt yr⁻¹ for a flux, or 0.0 if key is absent from budget."""
    da = budget.get(key)
    return 0.0 if da is None else _ann_Gt(_drop_placeholder_zero_members(s(da)))

def _fmt(v_Gt):
    v = abs(v_Gt)
    return f"{v/1e3:.1f}×10³" if v >= 500 else f"{v:.0f}"


# ── Budget QC: automated sanity checks ─────────────────────────────────────────
#
# Triage warnings for the same class of bugs documented in the README's
# Model-Specific Corrections table (wrong sign, unit/magnitude errors, per-member
# archiving bugs). Flags candidates for a human to check — never auto-corrects.

ICE_GAIN_TERMS  = ['basal growth', 'frazil', 'snowice']
ICE_LOSS_TERMS  = ['top melt', 'basal melt', 'lateral melt', 'evapsubl']
ICE_NET_TERMS   = ['dynamics']

SNOW_GAIN_TERMS = ['snowfall']
SNOW_LOSS_TERMS = ['snowmelt', 'evapsubl', 'snow to ice']
SNOW_NET_TERMS  = ['dynamics', 'wind drift']

def _term_annual_by_member(da, time_dim="time"):
    """Per-member annual-mean Gt yr⁻¹ (input: 10³ Gt month⁻¹), as a pandas Series
    indexed by member_id — no averaging across members (unlike _ann_Gt)."""
    if da is None or da.sizes.get(time_dim, 0) == 0 or da.sizes.get("member_id", 0) == 0:
        return None
    annual = da.groupby(f"{time_dim}.year").sum(time_dim)
    vals = (annual.mean("year") * 1e3).values
    return pd.Series(vals, index=annual["member_id"].values.astype(str))

def _member_model_prefix(member_ids):
    """Model source_id for each member_id string (everything before the first '_'),
    matching the convention used by sel_model/drop_sel."""
    return pd.Index(member_ids).str.replace(r'_.+', '', regex=True)

def _ratio_outliers(values, ratio_threshold):
    """Flag entries of `values` (a pandas Series of non-negative magnitudes) whose
    ratio to the group median exceeds ratio_threshold in either direction. Ratio-
    to-median rather than MAD/std since the known bugs here are large multiplicative
    errors (330x, 1800x, 1e6x) that would otherwise inflate their own spread stat."""
    flagged = {}
    if len(values) == 0:
        return flagged
    median = values.median()
    if median <= 0:
        return flagged
    for key, v in values.items():
        if v <= 0:
            continue
        ratio = max(v / median, median / v)
        if ratio > ratio_threshold:
            flagged[key] = (v, median, ratio)
    return flagged

def _qc_warn(warnings, verbose, region, budget_type, model, term, kind, detail):
    w = {"region": region, "budget": budget_type, "model": model, "term": term,
         "kind": kind, "detail": detail}
    warnings.append(w)
    if verbose:
        print(f"[QC] {region}/{budget_type} {model!r} {term!r}: {kind} — {detail}")
    return w

def run_budget_qc(budget, budget_type='ice', region_label='', outlier_ratio=5.0,
                   residual_frac=0.15, net_term_frac=0.10, min_group_size=3, verbose=True):
    """Triage a budget dict (as built in cmip6_sankey.ipynb, e.g. ice_budget_NH)
    for the categories of issue previously caught by manual inspection:

    - wrong_sign: a gain term averaging negative, or a loss term averaging positive,
      for some model (catches e.g. sidmassevapsubl reported loss-positive).
    - magnitude_outlier_model: a model's term magnitude is more than `outlier_ratio`x
      (or less than 1/`outlier_ratio`x) the cross-model median (catches unit/area-basis
      bugs like NorESM2's 330x snowfall).
    - magnitude_outlier_member: within an ensemble model (CESM2-WACCM, CESM2-LE, ...),
      one member's term magnitude is more than `outlier_ratio`x its sibling median
      (catches per-member archiving bugs like CESM2-WACCM's r1i1p1f1/r3i1p1f1 issues,
      which a model-mean check alone would miss).
    - closure_residual: gross gains vs. gross losses don't balance within
      `residual_frac` of total throughput (mirrors _balance_flows' Residual node).
    - net_term_nonzero: a term expected to net near-zero at full-region scale
      (dynamics, wind drift) exceeds `net_term_frac` of gross throughput.

    Returns the list of warning dicts (empty if nothing was flagged); also prints
    each warning when verbose=True.
    """
    gain_terms, loss_terms, net_terms = {
        'ice':  (ICE_GAIN_TERMS, ICE_LOSS_TERMS, ICE_NET_TERMS),
        'snow': (SNOW_GAIN_TERMS, SNOW_LOSS_TERMS, SNOW_NET_TERMS),
    }[budget_type]

    warnings = []

    per_member = {}
    for term, da in budget.items():
        vals = _term_annual_by_member(da)
        if vals is not None and len(vals) > 0:
            per_member[term] = vals

    if not per_member:
        return warnings

    any_index = next(iter(per_member.values())).index
    model_of = pd.Series(_member_model_prefix(any_index), index=any_index)

    model_means = {}  # term -> {model: mean_value}
    for term, vals in per_member.items():
        grouped = vals.groupby(model_of.reindex(vals.index))
        model_means[term] = grouped.mean()

        # 1. Sign convention, per model
        if term in gain_terms or term in loss_terms:
            expect_positive = term in gain_terms
            for model_name, mean_v in model_means[term].items():
                if expect_positive and mean_v < 0:
                    _qc_warn(warnings, verbose, region_label, budget_type, model_name, term,
                              "wrong_sign", f"expected positive (gain), got {mean_v:.0f} Gt/yr")
                elif not expect_positive and mean_v > 0:
                    _qc_warn(warnings, verbose, region_label, budget_type, model_name, term,
                              "wrong_sign", f"expected negative (loss), got {mean_v:.0f} Gt/yr")

        # 2a. Magnitude outliers across models
        abs_means = model_means[term].abs()
        if len(abs_means) >= min_group_size:
            for model_name, (v, median, ratio) in _ratio_outliers(abs_means, outlier_ratio).items():
                _qc_warn(warnings, verbose, region_label, budget_type, model_name, term,
                          "magnitude_outlier_model",
                          f"|{v:.0f}| Gt/yr is {ratio:.1f}x the cross-model median |{median:.0f}| Gt/yr")

        # 2b. Magnitude outliers across members within each ensemble model
        abs_vals = vals.abs()
        for model_name, member_ids in model_of.groupby(model_of).groups.items():
            group = abs_vals.reindex(member_ids).dropna()
            if len(group) < min_group_size:
                continue
            for member_id, (v, median, ratio) in _ratio_outliers(group, outlier_ratio).items():
                _qc_warn(warnings, verbose, region_label, budget_type, model_name, term,
                          "magnitude_outlier_member",
                          f"member {member_id!r}: |{v:.0f}| Gt/yr is {ratio:.1f}x the "
                          f"sibling median |{median:.0f}| Gt/yr")

    # 3 & 4. Closure residual and net-term-nonzero, per model
    all_models = sorted(model_of.unique())
    for model_name in all_models:
        gross_gain = sum(abs(model_means[t].get(model_name, 0.0)) for t in gain_terms if t in model_means)
        gross_loss = sum(abs(model_means[t].get(model_name, 0.0)) for t in loss_terms if t in model_means)
        total = gross_gain + gross_loss
        if total > 0:
            imbalance = gross_gain - gross_loss
            if abs(imbalance) / total > residual_frac:
                _qc_warn(warnings, verbose, region_label, budget_type, model_name, "(all terms)",
                          "closure_residual",
                          f"gains {gross_gain:.0f} vs. losses {gross_loss:.0f} Gt/yr "
                          f"({abs(imbalance)/total:.0%} of throughput)")
            for term in net_terms:
                if term not in model_means:
                    continue
                net_v = model_means[term].get(model_name, 0.0)
                if abs(net_v) / total > net_term_frac:
                    _qc_warn(warnings, verbose, region_label, budget_type, model_name, term,
                              "net_term_nonzero",
                              f"{net_v:.0f} Gt/yr is {abs(net_v)/total:.0%} of gross throughput "
                              f"(expected ≈0 at full-region scale)")

    return warnings


# ── Sankey: Plotly rendering ───────────────────────────────────────────────────

def _hex_to_rgba(hex_color, alpha=0.50):
    """Convert a '#rrggbb' hex color to a Plotly 'rgba(...)' string."""
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
                       hide_terminal_nodes=False, show_percent=False, show_group_percent=None,
                       sources_sinks_y=1.0, res_label_y=0.12, n_members=None):
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
    show_percent: if True, add each terminal node's share of its side's total
      (inflows as % of total sources, outflows as % of total sinks) in parentheses.
    show_group_percent: same, but for the combined Ocean Growth/Ocean Melt group
      node specifically. Defaults to show_percent — pass it explicitly when the
      terminal flow labels already carry their own custom text (e.g. the ±1 SD
      range in the uncertainty Sankeys) and shouldn't also get a plain % appended,
      but the group node — which has no such custom label — still should.
    subtitle: optional string shown as a small annotation in the top-left corner.
    """
    import plotly.graph_objects as go
    if show_group_percent is None:
        show_group_percent = show_percent

    inflows, outflows = _balance_flows(inflows, outflows)
    if not show_residual:
        inflows  = [(v, lb, col) for v, lb, col in inflows  if lb != "Residual"]
        outflows = [(v, lb, col) for v, lb, col in outflows if lb != "Residual"]
    n_in  = len(inflows)
    n_out = len(outflows)

    total_in  = sum(v for v, *_ in inflows)  or 1.0
    total_out = sum(v for v, *_ in outflows) or 1.0

    def _node_label(lb, v, total=None, want_percent=None):
        want_percent = show_percent if want_percent is None else want_percent
        parts = []
        if show_values:
            parts.append(f"{_fmt(v)} Gt yr⁻¹")
        if want_percent and total is not None:
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

    def _in_group(lb, members):
        """True if lb is (or starts with, for percentage-annotated uncertainty
        labels like "Basal growth 45% (39-51%)") one of the group's member terms."""
        return any(lb == m or lb.startswith(m + " ") for m in members)

    group_idx = None
    if inflow_group is not None:
        members = inflow_group["members"]
        grouped = [i for i, (*_, lb, _) in enumerate(inflows) if _in_group(lb, members)]
        if grouped:
            combined_total = sum(inflows[i][0] for i in grouped)
            combined_y = sum(node_y[i] * inflows[i][0] for i in grouped) / combined_total
            group_idx = len(node_labels)
            group_name = inflow_group.get("name")
            group_label = (_node_label(group_name or "", combined_total, total_in, want_percent=show_group_percent)
                           if (show_values or show_group_percent) else (group_name or ""))
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
        out_grouped = [j for j, (*_, lb, _) in enumerate(outflows) if _in_group(lb, out_members)]
        if out_grouped:
            combined_out_total = sum(outflows[j][0] for j in out_grouped)
            combined_out_y = sum(out_ys[j] * outflows[j][0] for j in out_grouped) / combined_out_total
            out_group_idx = len(node_labels)
            out_group_name = outflow_group.get("name")
            out_group_label = (_node_label(out_group_name or "", combined_out_total, total_out, want_percent=show_group_percent)
                               if (show_values or show_group_percent) else (out_group_name or ""))
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
        target = group_idx if (group_idx is not None and _in_group(lb, inflow_group["members"])) else res_idx
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
        source = out_group_idx if (out_group_idx is not None and _in_group(lb, outflow_group["members"])) else res_idx
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

    # The reservoir node's own label is pulled out of the inline Sankey node label (which
    # Plotly draws right next to the node, crowding the adjacent source/sink flow labels)
    # and rendered as a standalone annotation above the plot instead. Hover text is kept
    # intact via customdata. The "Ocean Melt"/"Ocean Growth" grouping node gets the same
    # treatment, but anchored to the left of the node instead, so it doesn't collide with
    # the terminal sink labels to its right.
    node_display_labels = list(node_labels)
    node_display_labels[res_idx] = ""
    if out_group_idx is not None:
        node_display_labels[out_group_idx] = ""

    NODE_THICKNESS = 14
    NODE_PAD = 22
    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=NODE_PAD,
            thickness=NODE_THICKNESS,
            line=dict(color="white", width=0),
            label=node_display_labels,
            color=node_colors,
            x=node_x,
            y=node_y,
            customdata=node_labels,
            hovertemplate="%{customdata}<extra></extra>",
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
    fig_id = f"{budget_type} ({model_label}{', ' + subtitle if subtitle else ''})"
    if scale_ref is not None and scale_ref[0] > 0:
        ref_total, ref_n_busiest = scale_ref
        reservoir_total = sum(v for v, *_ in inflows)
        gaps_this = NODE_PAD * max(0, n_busiest - 1)
        gaps_ref  = NODE_PAD * max(0, ref_n_busiest - 1)
        available = (PLOT_AREA_FULL - gaps_ref) * (reservoir_total / ref_total)
        plot_area = max(50, gaps_this + available)
        print(f"[normalization] {fig_id}: reservoir={reservoir_total:,.1f} Gt yr⁻¹ / "
              f"scale_ref={ref_total:,.1f} Gt yr⁻¹ -> scale factor={reservoir_total / ref_total:.3f}, "
              f"plot_area={plot_area:.0f}px")
    else:
        plot_area = PLOT_AREA_FULL
        print(f"[normalization] {fig_id}: no scale_ref given, using full plot area ({plot_area}px)")
    height = MARGIN_T + MARGIN_B + plot_area

    if title is not None:
        title_text = title
    else:
        model_id = "MME" if model_label == "Multi-model mean" else model_label
        if n_members is not None and model_id != "MME":
            member_word = "member" if n_members == 1 else "members"
            model_id = f"{model_id}, {n_members} {member_word}"
        prefix = f"{subtitle} " if subtitle else ""
        title_text = f"{prefix}{budget_type} Budget ({model_id})"

    annotations = [
        dict(x=node_x[res_idx], y=res_label_y, xref="paper", yref="paper",
             text=res_label, showarrow=False,
             xanchor="center", yanchor="top",
             font=dict(size=11, color="#333333")),
        dict(x=0.47, y=sources_sinks_y, xref="paper", yref="paper", text="sources",
             showarrow=False, xanchor="right",
             font=dict(size=11, color="#888888")),
        dict(x=0.53, y=sources_sinks_y, xref="paper", yref="paper", text="sinks",
             showarrow=False, xanchor="left",
             font=dict(size=11, color="#888888")),
    ]
    if out_group_idx is not None:
        # node.y (domain space, increases downward) and this annotation's y (paper
        # space, increases upward) map to the same plot-area pixel range but run in
        # opposite directions, so tracking the node vertically is a straight flip.
        annotations.append(
            dict(x=node_x[out_group_idx], y=1 - node_y[out_group_idx], xref="paper", yref="paper",
                 text=node_labels[out_group_idx], showarrow=False,
                 xanchor="right", yanchor="middle", xshift=-8,
                 font=dict(size=11, color="#333333"))
        )

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
        annotations=annotations,
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
        (abs(si_v), "Snow-to-ice",    _C["snow2ice_g"]),
        (abs(fr_v), "Open water ice production",      _C["frazil"]),
        (abs(bg_v), "Basal growth", _C["basal_g"]),
        (dyn_in,    "Dynamics",  _C["dyn_in"]),
    ]
    outflows = [
        (abs(es_v), "Vapour Exchange",   _C["evapsubl"]),
        (abs(tm_v), "Top melt",     _C["top_melt"]),
        (abs(lm_v), "Lateral melt", _C["lat_melt"]),
        (abs(bm_v), "Basal melt",   _C["basal_melt"]),
        (dyn_out,   "Dynamics",  _C["dyn_out"]),
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
    n_members = s(ice_budget["basal growth"]).sizes["member_id"] if model is not None else None
    if model is not None and n_members == 0:
        return go.Figure()

    inflows, outflows = _ice_flows(ice_budget, model)
    inflow_group  = {"members": {"Open water ice production", "Basal growth"}, "color": _C["ice_growth"], "name": inflow_group_name}
    outflow_group = {"members": {"Basal melt", "Lateral melt"}, "color": _C["ice_melt"], "name": outflow_group_name}
    label = model if model is not None else "Multi-model mean"
    return _plotly_sankey_fig(res_label or "Sea Ice", res_color or _C["ice_res"],
                              inflows, outflows, label, "Sea Ice Mass", title=title,
                              inflow_group=inflow_group, outflow_group=outflow_group, scale_ref=scale_ref,
                              show_values=show_values, show_residual=show_residual, subtitle=subtitle,
                              hide_terminal_nodes=hide_terminal_nodes, show_percent=show_percent,
                              n_members=n_members)


# ── Sankey: uncertainty (inter-ensemble-member spread) ─────────────────────────
#
# Adds a "how uncertain is this flow" layer on top of the existing MME Sankeys.
# Spread is the inter-member standard deviation (Gt yr⁻¹) pooled across the
# ensemble members of every model in the budget dict — ensemble-member spread,
# not (yet) decomposed into inter-model spread. Reported as the range a term's
# own % contribution to the budget would span if the term alone were shifted by
# ±1 SD. Rendering reuses _plotly_sankey_fig, so styling matches the plain Sankeys.

def _ann_Gt_per_member(da, member_dim="member_id", time_dim="time"):
    """Per-member annual Gt yr⁻¹ (input: 10³ Gt month⁻¹) — one value per member_id,
    unlike _ann_Gt which averages across members."""
    annual = da.groupby(f"{time_dim}.year").sum(time_dim)
    return annual.mean("year") * 1e3

def _term_spread(budget, key, model=None):
    """Mean and standard deviation (both Gt yr⁻¹) for one budget term, computed from
    per-member annual values pooled across every model in `budget` (or restricted to
    one model's members if `model` is given)."""
    da = budget.get(key)
    if da is None:
        return 0.0, 0.0
    da_sel = _drop_placeholder_zero_members(_sel(da, model))
    if da_sel.sizes.get("member_id", 0) < 2 or da_sel.sizes.get("time", 0) == 0:
        return _ann_Gt(da_sel), 0.0
    per_member = _ann_Gt_per_member(da_sel)
    mean = float(per_member.mean("member_id").values)
    std = float(per_member.std("member_id").values)
    return mean, std

def _ice_flows_uncertainty(ice_budget, model=None):
    """Like _ice_flows, but each (value, label, color) tuple gains a fourth element:
    the term's inter-member standard deviation in Gt yr⁻¹ (see _term_spread)."""
    terms = {
        "basal growth": ("Basal growth", _C["basal_g"]),
        "frazil":       ("Open water ice production", _C["frazil"]),
        "snowice":      ("Snow-to-ice", _C["snow2ice_g"]),
        "top melt":     ("Top melt", _C["top_melt"]),
        "basal melt":   ("Basal melt", _C["basal_melt"]),
        "lateral melt": ("Lateral melt", _C["lat_melt"]),
        "evapsubl":     ("Vapour Exchange", _C["evapsubl"]),
    }
    v = {key: (*_term_spread(ice_budget, key, model), label, color)
         for key, (label, color) in terms.items()}
    # v[key] = (mean, std, label, color)

    dy_mean, dy_std = _term_spread(ice_budget, "dynamics", model)

    inflows = [
        (abs(v["snowice"][0]),      v["snowice"][2],      v["snowice"][3],      v["snowice"][1]),
        (abs(v["frazil"][0]),       v["frazil"][2],       v["frazil"][3],       v["frazil"][1]),
        (abs(v["basal growth"][0]), v["basal growth"][2], v["basal growth"][3], v["basal growth"][1]),
        (max(0.0, dy_mean),         "Dynamics",           _C["dyn_in"],         dy_std),
    ]
    outflows = [
        (abs(v["evapsubl"][0]),     v["evapsubl"][2],     v["evapsubl"][3],     v["evapsubl"][1]),
        (abs(v["top melt"][0]),     v["top melt"][2],     v["top melt"][3],     v["top melt"][1]),
        (abs(v["lateral melt"][0]), v["lateral melt"][2], v["lateral melt"][3], v["lateral melt"][1]),
        (abs(v["basal melt"][0]),   v["basal melt"][2],   v["basal melt"][3],   v["basal melt"][1]),
        (max(0.0, -dy_mean),        "Dynamics",           _C["dyn_out"],        dy_std),
    ]
    return inflows, outflows

def _hex_blend_white(hex_color, frac):
    """Blend hex_color toward white by frac (0 = full color, 1 = white)."""
    h = hex_color.lstrip('#')
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    r, g, b = (int(c + (255 - c) * frac) for c in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"

def _ice_sankey_uncertainty_flows(ice_budget, model=None, include_dynamics=True,
                                  shade_uncertainty=False, max_rel_spread=0.6):
    """Build (inflows, outflows) as (value, label, color) 3-tuples ready for
    _plotly_sankey_fig: `value` is each term's ensemble-mean Gt yr⁻¹, `label` shows
    the term's mean % share of its side's total plus the range that share would span
    if the term alone were shifted by ±1 SD (holding the total and every other term
    fixed at their means) — i.e. how much the ensemble's disagreement on this one term
    could move its own contribution.

    `color` is the term's plain base color by default. Set shade_uncertainty=True to
    also fade it toward white in proportion to relative spread (std/mean, capped at
    max_rel_spread) — off by default since it reads as harder to parse at a glance
    than the label's ±1 SD range alone.

    include_dynamics: set False to drop the Dynamics flow entirely (and exclude it
    from the % totals) — for a full hemisphere (SH/NH) this term is a near-zero
    residual of the (closed) integration domain rather than a genuine export, so it's
    dropped there; the Inner Arctic (IA) domain is where dynamics is a real signal.
    """
    inflows_raw, outflows_raw = _ice_flows_uncertainty(ice_budget, model)
    inflows_raw  = [(v, lb, col, sd) for v, lb, col, sd in inflows_raw  if v > 0]
    outflows_raw = [(v, lb, col, sd) for v, lb, col, sd in outflows_raw if v > 0]
    if not include_dynamics:
        inflows_raw  = [t for t in inflows_raw  if t[1] != "Dynamics"]
        outflows_raw = [t for t in outflows_raw if t[1] != "Dynamics"]
    total_in  = sum(v for v, *_ in inflows_raw)  or 1.0
    total_out = sum(v for v, *_ in outflows_raw) or 1.0

    def _build(terms, total):
        built = []
        for v, lb, col, sd in terms:
            pct    = 100 * v / total
            pct_lo = 100 * max(0.0, v - sd) / total
            pct_hi = 100 * (v + sd) / total
            if shade_uncertainty:
                rel = min(sd / v, max_rel_spread) / max_rel_spread if v > 0 else 0.0
                col = _hex_blend_white(col, 0.6 * rel)
            label = f"{lb} {pct:.0f}% ({pct_lo:.0f}–{pct_hi:.0f}%)"
            built.append((v, label, col))
        return built

    return _build(inflows_raw, total_in), _build(outflows_raw, total_out)

def make_ice_sankey_uncertainty_plotly(ice_budget, model=None, res_label="Sea Ice",
                                       res_color=None, title=None, subtitle=None,
                                       inflow_group_name=None, outflow_group_name=None,
                                       scale_ref=None, show_residual=True, width=1100,
                                       include_dynamics=True,
                                       shade_uncertainty=False, max_rel_spread=0.6):
    """Ice mass budget Sankey with its node labels showing the ±1 SD range of each
    flow's % contribution to the budget (see _ice_sankey_uncertainty_flows). Rendered
    via _plotly_sankey_fig with hide_terminal_nodes=True, so it matches
    make_ice_sankey_plotly's styling (no colored source/sink bars, reservoir label as
    an annotation, same margins, same Ocean Growth / Ocean Melt grouping node when
    inflow_group_name / outflow_group_name are given) — only the flow labels differ.

    include_dynamics: see _ice_sankey_uncertainty_flows — pass False for full
    hemisphere (SH/NH) budgets, True for the Inner Arctic (IA).
    shade_uncertainty / max_rel_spread: set shade_uncertainty=True to also fade each
    flow's ribbon color toward white by relative spread (off by default — this reads
    as harder to parse at a glance than the label's ±1 SD range alone).
    """
    inflows, outflows = _ice_sankey_uncertainty_flows(ice_budget, model=model,
                                                       include_dynamics=include_dynamics,
                                                       shade_uncertainty=shade_uncertainty,
                                                       max_rel_spread=max_rel_spread)
    inflow_group  = {"members": {"Open water ice production", "Basal growth"}, "color": _C["ice_growth"], "name": inflow_group_name}
    outflow_group = {"members": {"Basal melt", "Lateral melt"}, "color": _C["ice_melt"], "name": outflow_group_name}
    label = model if model is not None else "Multi-model mean"
    n_members = _sel(ice_budget["basal growth"], model).sizes["member_id"] if model is not None else None
    fig = _plotly_sankey_fig(res_label, res_color or _C["ice_res"], inflows, outflows,
                             label, "Sea Ice Mass", title=title,
                             inflow_group=inflow_group, outflow_group=outflow_group, scale_ref=scale_ref,
                             show_values=False, show_residual=show_residual,
                             subtitle=subtitle, hide_terminal_nodes=True,
                             show_percent=False, show_group_percent=True, n_members=n_members)
    fig.update_layout(width=width)
    return fig


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
        (dyn_in,    "Dynamics", _C["dyn_in"]),
        (wind_in,   "Wind drift", _C["wind_in"]),
    ]
    outflows = [
        (abs(es_v),  "Vapour Exchange",  _C["evapsubl"]),
        (wind_out,   "Wind drift", _C["wind_out"]),
        (abs(sm_v),  "Snowmelt",    _C["snowmelt"]),
        (abs(s2_v),  "Snow-to-ice",   _C["snow2ice_s"]),
        (dyn_out,    "Dynamics", _C["dyn_out"]),
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
    n_members = s(snow_budget["snowfall"]).sizes["member_id"] if model is not None else None
    if model is not None and n_members == 0:
        return go.Figure()

    inflows, outflows = _snow_flows(snow_budget, model)
    label = model if model is not None else "Multi-model mean"
    return _plotly_sankey_fig(res_label or "Snow", res_color or _C["snow_res"],
                              inflows, outflows, label, "Snow on Sea Ice Mass", title=title,
                              scale_ref=scale_ref, show_values=show_values,
                              show_residual=show_residual, subtitle=subtitle,
                              hide_terminal_nodes=hide_terminal_nodes, show_percent=show_percent,
                              sources_sinks_y=1.03, res_label_y=0.06, n_members=n_members)


# ── Sankey: snow uncertainty (inter-ensemble-member spread) ────────────────────
#
# Mirrors the ice uncertainty section above (see its comment for the Color-code
# technique this follows) — same ribbon-fading-by-spread and ±1 SD node-label
# treatment, applied to the snow budget's terms instead.

def _snow_flows_uncertainty(snow_budget, model=None):
    """Like _snow_flows, but each (value, label, color) tuple gains a fourth element:
    the term's inter-member standard deviation in Gt yr⁻¹ (see _term_spread)."""
    terms = {
        "snowfall":    ("Snowfall", _C["snowfall"]),
        "snowmelt":    ("Snowmelt", _C["snowmelt"]),
        "snow to ice": ("Snow-to-ice", _C["snow2ice_s"]),
        "evapsubl":    ("Vapour Exchange", _C["evapsubl"]),
    }
    v = {key: (*_term_spread(snow_budget, key, model), label, color)
         for key, (label, color) in terms.items()}
    # v[key] = (mean, std, label, color)

    dy_mean, dy_std = _term_spread(snow_budget, "dynamics", model)
    wd_mean, wd_std = _term_spread(snow_budget, "wind drift", model)

    inflows = [
        (abs(v["snowfall"][0]), v["snowfall"][2], v["snowfall"][3], v["snowfall"][1]),
        (max(0.0, dy_mean),     "Dynamics",       _C["dyn_in"],     dy_std),
        (max(0.0, wd_mean),     "Wind drift",     _C["wind_in"],    wd_std),
    ]
    outflows = [
        (abs(v["evapsubl"][0]),    v["evapsubl"][2],    v["evapsubl"][3],    v["evapsubl"][1]),
        (max(0.0, -wd_mean),       "Wind drift",        _C["wind_out"],      wd_std),
        (abs(v["snowmelt"][0]),    v["snowmelt"][2],    v["snowmelt"][3],    v["snowmelt"][1]),
        (abs(v["snow to ice"][0]), v["snow to ice"][2], v["snow to ice"][3], v["snow to ice"][1]),
        (max(0.0, -dy_mean),       "Dynamics",          _C["dyn_out"],       dy_std),
    ]
    return inflows, outflows

def _snow_sankey_uncertainty_flows(snow_budget, model=None, shade_uncertainty=False, max_rel_spread=0.6):
    """Build (inflows, outflows) as (value, label, color) 3-tuples ready for
    _plotly_sankey_fig — see _ice_sankey_uncertainty_flows for what `value`/`color`/
    `label`/`shade_uncertainty`/`max_rel_spread` mean here. Terms with zero mean (e.g.
    wind drift/dynamics, archived by only a handful of models) are dropped by the same
    v > 0 filter as the plain make_snow_sankey_plotly, so no separate
    include_dynamics/include_wind_drift toggle is needed here."""
    inflows_raw, outflows_raw = _snow_flows_uncertainty(snow_budget, model)
    inflows_raw  = [(v, lb, col, sd) for v, lb, col, sd in inflows_raw  if v > 0]
    outflows_raw = [(v, lb, col, sd) for v, lb, col, sd in outflows_raw if v > 0]
    total_in  = sum(v for v, *_ in inflows_raw)  or 1.0
    total_out = sum(v for v, *_ in outflows_raw) or 1.0

    def _build(terms, total):
        built = []
        for v, lb, col, sd in terms:
            pct    = 100 * v / total
            pct_lo = 100 * max(0.0, v - sd) / total
            pct_hi = 100 * (v + sd) / total
            if shade_uncertainty:
                rel = min(sd / v, max_rel_spread) / max_rel_spread if v > 0 else 0.0
                col = _hex_blend_white(col, 0.6 * rel)
            label = f"{lb} {pct:.0f}% ({pct_lo:.0f}–{pct_hi:.0f}%)"
            built.append((v, label, col))
        return built

    return _build(inflows_raw, total_in), _build(outflows_raw, total_out)

def make_snow_sankey_uncertainty_plotly(snow_budget, model=None, res_label="Snow",
                                        res_color=None, title=None, subtitle=None,
                                        scale_ref=None, show_residual=True, width=1100,
                                        shade_uncertainty=False, max_rel_spread=0.6):
    """Snow mass budget Sankey with its node labels showing the ±1 SD range of each
    flow's % contribution to the budget (see _snow_sankey_uncertainty_flows). Rendered
    via _plotly_sankey_fig to match make_snow_sankey_plotly's styling — only the flow
    labels differ.

    shade_uncertainty / max_rel_spread: see make_ice_sankey_uncertainty_plotly."""
    inflows, outflows = _snow_sankey_uncertainty_flows(snow_budget, model=model,
                                                        shade_uncertainty=shade_uncertainty,
                                                        max_rel_spread=max_rel_spread)
    label = model if model is not None else "Multi-model mean"
    n_members = _sel(snow_budget["snowfall"], model).sizes["member_id"] if model is not None else None
    fig = _plotly_sankey_fig(res_label, res_color or _C["snow_res"], inflows, outflows,
                             label, "Snow on Sea Ice Mass", title=title,
                             scale_ref=scale_ref, show_values=False, show_residual=show_residual,
                             subtitle=subtitle, hide_terminal_nodes=True, show_percent=False,
                             sources_sinks_y=1.03, res_label_y=0.06, n_members=n_members)
    fig.update_layout(width=width)
    return fig


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
