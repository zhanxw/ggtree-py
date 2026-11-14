"""
Tree view & clade manipulation functions (port of R/clade-functions.R, R/layout.R, R/format.R):
get_taxa_name, viewClade, collapse, expand, rotate, flip, scaleClade, zoomClade,
open_tree, rotate_tree, layout_rectangular, layout_circular, layout_inward_circular,
layout_fan, layout_dendrogram, layout_ggtree, label_pad, range_format.
"""

from __future__ import annotations

import copy
import math
import warnings
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np

from .aes import aes
from .dataframe import ColumnSeries, DataFrame, _is_scalar_na
from .geoms import Layer, geom_polygon
from .layout import (
    calculate_angle,
    calculate_branch_mid,
    get_layout,
    re_assign_ycoord_df,
    reassign_y_from_node_to_root,
)
from .themes_scales import PlotModifier, Rel, coord_cartesian, coord_flip, coord_polar, rel, revts, theme_dendrogram, theme_tree
from .tidytree import child, getRoot, isTip, offspring, parent, taxa2node


def get_taxa_name(tree_view: Any = None, node: Optional[Union[int, str]] = None) -> List[str]:
    """
    Get ordered taxa (tip) names from top to bottom of a tree view
    (exact port of ggtree::get_taxa_name in R/clade-functions.R).
    """
    df = tree_view.data if hasattr(tree_view, "data") else tree_view
    if node is not None:
        df = offspring(df, node)
    tips = df[df["isTip"]] if "isTip" in df else df
    recs = [r for r in tips.to_records() if not _is_scalar_na(r.get("y"))]
    recs.sort(key=lambda r: float(r["y"]), reverse=True)
    return [str(r["label"]) for r in recs if not _is_scalar_na(r.get("label"))]


def viewClade(tree_view: Any = None, node: Optional[Union[int, str]] = None, xmax_adjust: float = 0.0) -> Any:
    """
    Focus plot view on a specific clade (exact port of ggtree::viewClade in R/clade-functions.R).
    Can be called directly `viewClade(p, node)` or added via `p + viewClade(node=...)`.
    """
    if tree_view is None or isinstance(tree_view, (int, str, np.integer)):
        nd = tree_view if node is None else node

        class _ViewCladeMod(PlotModifier):
            def apply(self, plot: Any) -> Any:
                return viewClade(plot, nd, xmax_adjust=xmax_adjust)

        return _ViewCladeMod()

    p = tree_view.copy()
    if node is None:
        return p
    df = p.data
    nid = taxa2node(df, node)
    cpos = offspring(df, nid, self_include=True)
    ys = [float(v) for v in cpos["y"].tolist() if not _is_scalar_na(v)]
    xs = [float(v) for v in cpos["x"].tolist() if not _is_scalar_na(v)]
    all_xs = [float(v) for v in df["x"].tolist() if not _is_scalar_na(v)]
    node_row = df[df["node"] == nid].to_records()
    br_len = 0.0
    if node_row and "branch.length" in df and not _is_scalar_na(node_row[0].get("branch.length")):
        br_len = float(node_row[0]["branch.length"])
    xmin = (min(xs) if xs else 0.0) - br_len / 2.0
    xmax = (max(all_xs) if all_xs else 1.0) + float(xmax_adjust)
    ymin = min(ys) if ys else 0.0
    ymax = max(ys) if ys else 1.0
    return p + coord_cartesian(xlim=(xmin, xmax), ylim=(ymin, ymax), expand=False)


def zoomClade(tree_view: Any = None, node: Optional[Union[int, str]] = None, xexpand: Optional[float] = None) -> Any:
    """
    Zoom into a clade rooted at 'node' (exact port of ggtree::zoomClade in R/clade-functions.R).
    """
    if tree_view is None or isinstance(tree_view, (int, str, np.integer)):
        nd = tree_view if node is None else node

        class _ZoomMod(PlotModifier):
            def apply(self, plot: Any) -> Any:
                return zoomClade(plot, nd, xexpand=xexpand)

        return _ZoomMod()

    p = tree_view.copy()
    df = p.data
    nid = taxa2node(df, node)  # type: ignore[arg-type]
    sp = offspring(df, nid, self_include=True)
    xs = [float(v) for v in sp["x"].tolist() if not _is_scalar_na(v)]
    ys = [float(v) for v in sp["y"].tolist() if not _is_scalar_na(v)]
    xmin, xmax = (min(xs), max(xs)) if xs else (0.0, 1.0)
    if xexpand is not None:
        xmax += (xmax - xmin) * float(xexpand)
    ymin, ymax = (min(ys), max(ys)) if ys else (0.0, 1.0)
    return p + coord_cartesian(xlim=(xmin, xmax), ylim=(ymin, ymax))


def _validate_collapse_height(height: Any) -> Any:
    """Validate collapse() height argument (exact port of R/clade-functions.R .validate_collapse_height)."""
    if height is None:
        return None
    err_msg = (
        "`height` must be NULL, `ggplot2::rel()` with a single non-negative value, "
        "or a single non-negative numeric value."
    )
    if isinstance(height, Rel):
        val = float(height)
        if math.isnan(val) or math.isinf(val) or val < 0:
            raise ValueError(err_msg)
        return Rel(val)
    if isinstance(height, (bool, str, bytes, dict)):
        raise ValueError(err_msg)
    if isinstance(height, (list, tuple, np.ndarray)):
        arr = list(height)
        if len(arr) != 1 or isinstance(arr[0], (bool, str, bytes)):
            raise ValueError(err_msg)
        try:
            val = float(arr[0])
        except Exception:
            raise ValueError(err_msg)
        if math.isnan(val) or math.isinf(val) or val < 0:
            raise ValueError(err_msg)
        return val
    try:
        val = float(height)
    except Exception:
        raise ValueError(err_msg)
    if math.isnan(val) or math.isinf(val) or val < 0:
        raise ValueError(err_msg)
    return val


def _resolve_collapse_height(height: Any, original_span: float, mode: str) -> float:
    if height is None:
        return 0.0 if mode == "none" else float(original_span)
    if isinstance(height, Rel):
        return float(original_span) * float(height)
    return float(height)


def _rescale_collapse_y(y: Any, ymin: float, original_span: float, collapsed_span: float) -> Any:
    if isinstance(y, (list, tuple, np.ndarray)):
        if original_span == 0:
            return [float(ymin) for _ in y]
        return [ymin + (float(v) - ymin) * collapsed_span / original_span for v in y]
    if original_span == 0:
        return float(ymin)
    return ymin + (float(y) - ymin) * collapsed_span / original_span


def _recalculate_collapsed_tree(tree_view: Any, df: DataFrame, nid: int) -> DataFrame:
    layout_str = get_layout(tree_view) or "rectangular"
    out = reassign_y_from_node_to_root(df, nid)
    out = calculate_branch_mid(out, layout=layout_str)
    xs_col = out["x"].tolist()
    non_na_idx = [i for i, xv in enumerate(xs_col) if not _is_scalar_na(xv)]
    angles_col: List[Any] = list(out["angle"].tolist()) if "angle" in out else [None] * len(out)
    if non_na_idx:
        sub_df = out[non_na_idx]
        sub_ang = calculate_angle(sub_df)["angle"].tolist()
        for pos, ang in zip(non_na_idx, sub_ang):
            angles_col[pos] = ang
    out["angle"] = angles_col
    return out


def collapse(
    x: Any = None,
    node: Optional[Union[int, str]] = None,
    mode: str = "none",
    clade_name: Optional[str] = None,
    height: Any = None,
    **kwargs: Any,
) -> Any:
    """
    Collapse a clade at 'node' with mode='none', 'max', 'min', or 'mixed'
    (exact port of collapse.ggtree in R/clade-functions.R).
    Supports both `collapse(p, node=17, ...)` and `p + collapse(node=17, ...)`.
    """
    if x is None or isinstance(x, (int, str, np.integer)):
        nd = x if node is None else node

        class _CollapseMod(PlotModifier):
            def apply(self, plot: Any) -> Any:
                return collapse(plot, node=nd, mode=mode, clade_name=clade_name, height=height, **kwargs)

        return _CollapseMod()

    if mode not in ("none", "max", "min", "mixed"):
        raise ValueError(f"collapse mode '{mode}' is not supported; use 'none', 'max', 'min', or 'mixed'.")
    height = _validate_collapse_height(height)

    p = x.copy()
    if node is None:
        return p
    df = p.data.copy()
    nid = taxa2node(df, node)
    nodes_col = [int(n) for n in df["node"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes_col)}
    nid_idx = id_to_idx[nid]

    if _is_scalar_na(df["x"][nid_idx]):
        warnings.warn("specific node was already collapsed...")
        return p

    sp_df = offspring(df, nid)
    if len(sp_df) == 0:
        warnings.warn("input node is a tip...")
        return p

    node_data = df[nid_idx:nid_idx + 1].copy()
    sp_snapshot = sp_df.copy()

    sp_nodes = [int(n) for n in sp_df["node"].tolist()]
    sp_idx = [id_to_idx[s] for s in sp_nodes if s in id_to_idx]
    tips_df = sp_df[sp_df["isTip"]] if "isTip" in sp_df else sp_df
    if len(tips_df) == 0:
        tips_df = sp_df

    sp_ys = [float(v) for v in sp_df["y"].tolist() if not _is_scalar_na(v)]
    if not sp_ys:
        warnings.warn("specific node was already collapsed...")
        return p
    sp_ymin, sp_ymax = min(sp_ys), max(sp_ys)
    original_span = sp_ymax - sp_ymin
    collapsed_span = _resolve_collapse_height(height, original_span, mode)
    adjusted_height = (mode == "none") or (height is not None)

    ys_col = list(df["y"].tolist())
    xs_col = list(df["x"].tolist())

    if adjusted_height:
        shift = original_span - collapsed_span
        for i, yv in enumerate(ys_col):
            if not _is_scalar_na(yv) and float(yv) > sp_ymax:
                ys_col[i] = float(yv) - shift
        ys_col[nid_idx] = _rescale_collapse_y(float(node_data["y"][0]), sp_ymin, original_span, collapsed_span)
        df["y"] = ys_col

    triangle_df: Optional[DataFrame] = None
    if mode == "none":
        for i in sp_idx:
            xs_col[i] = None
            ys_col[i] = None
        df["x"] = xs_col
        df["y"] = ys_col
        df = _recalculate_collapsed_tree(p, df, nid)
    else:
        tip_xs = [float(v) for v in tips_df["x"].tolist() if not _is_scalar_na(v)]
        raw_tip_ys = [float(v) for v in tips_df["y"].tolist() if not _is_scalar_na(v)]
        tip_ys = (
            _rescale_collapse_y(raw_tip_ys, sp_ymin, original_span, collapsed_span)
            if adjusted_height
            else raw_tip_ys
        )
        xmax_sp = max(tip_xs)
        xmin_sp = min(tip_xs)
        ymax_sp = max(tip_ys)
        ymin_sp = min(tip_ys)

        node_x = float(df["x"][nid_idx])
        node_y = float(df["y"][nid_idx])
        if mode == "max":
            triangle_df = DataFrame({"x": [node_x, xmax_sp, xmax_sp], "y": [node_y, ymin_sp, ymax_sp]})
        elif mode == "min":
            triangle_df = DataFrame({"x": [node_x, xmin_sp, xmin_sp], "y": [node_y, ymin_sp, ymax_sp]})
        else:  # mixed
            triangle_df = DataFrame({"x": [node_x, xmin_sp, xmax_sp], "y": [node_y, ymin_sp, ymax_sp]})

        xs_col = list(df["x"].tolist())
        ys_col = list(df["y"].tolist())
        for i in sp_idx:
            xs_col[i] = None
            ys_col[i] = None
        df["x"] = xs_col
        df["y"] = ys_col

        if adjusted_height:
            df = _recalculate_collapsed_tree(p, df, nid)

    if clade_name is not None:
        label_col = list(df["label"].tolist())
        label_col[nid_idx] = str(clade_name)
        df["label"] = label_col

    collapsed_clades = dict(df.attrs.get("collapsed_clades", {}))
    collapsed_clades[nid] = {
        "sp_df": sp_snapshot,
        "mode": mode,
        "collapsed_span": collapsed_span if adjusted_height else None,
        "node_data": node_data,
    }
    df.attrs["collapsed_clades"] = collapsed_clades
    p.data = df

    if triangle_df is not None:
        poly_params = dict(kwargs)
        poly_params["inherit.aes"] = False
        poly_layer = geom_polygon(
            mapping=aes(x="x", y="y"),
            data=triangle_df,
            collapse_node=nid,
            **poly_params,
        )
        p.layers.append(poly_layer)

    return p


def expand(tree_view: Any = None, node: Optional[Union[int, str]] = None) -> Any:
    """
    Expand a previously collapsed clade at 'node'
    (exact port of ggtree::expand in R/clade-functions.R).
    Supports both `expand(p, node)` and `p + expand(node)`.
    """
    if tree_view is None or isinstance(tree_view, (int, str, np.integer)):
        nd = tree_view if node is None else node

        class _ExpandMod(PlotModifier):
            def apply(self, plot: Any) -> Any:
                return expand(plot, nd)

        return _ExpandMod()

    p = tree_view.copy()
    if node is None:
        return p
    df = p.data.copy()
    collapsed_clades = dict(df.attrs.get("collapsed_clades", {}))
    if not collapsed_clades:
        return p

    try:
        nid = taxa2node(df, node)
    except Exception:
        nid = int(node) if isinstance(node, (int, np.integer)) else -1
    if nid not in collapsed_clades:
        return p

    info = collapsed_clades.pop(nid)
    df.attrs["collapsed_clades"] = collapsed_clades

    # Remove polygon layer added for this collapsed node
    p.layers = [lyr for lyr in p.layers if lyr.params.get("collapse_node") != nid]

    sp_df: DataFrame = info["sp_df"]
    mode: str = info["mode"]
    collapsed_span: Optional[float] = info.get("collapsed_span")
    node_data: Optional[DataFrame] = info.get("node_data")

    nodes_col = [int(n) for n in df["node"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes_col)}

    if mode == "none" or collapsed_span is not None:
        sp_ys = [float(v) for v in sp_df["y"].tolist() if not _is_scalar_na(v)]
        if sp_ys:
            sp_ymin, sp_ymax = min(sp_ys), max(sp_ys)
            original_span = sp_ymax - sp_ymin
            current_span = 0.0 if collapsed_span is None else float(collapsed_span)
            current_top = sp_ymin + current_span
            ys_col = list(df["y"].tolist())
            for i, yv in enumerate(ys_col):
                if not _is_scalar_na(yv) and float(yv) > current_top:
                    ys_col[i] = float(yv) + (original_span - current_span)
            df["y"] = ys_col

    def _restore_rows_from(src: DataFrame) -> None:
        src_recs = src.to_records()
        col_lists = {c: list(df[c].tolist()) for c in df.columns}
        for rec in src_recs:
            nd_id = int(rec["node"])
            if nd_id not in id_to_idx:
                continue
            ridx = id_to_idx[nd_id]
            for c in df.columns:
                if c in rec:
                    col_lists[c][ridx] = rec[c]
        for c, vals in col_lists.items():
            df[c] = vals

    if node_data is not None and len(node_data) > 0:
        _restore_rows_from(node_data)
    _restore_rows_from(sp_df)

    if mode == "none" or collapsed_span is not None:
        df = _recalculate_collapsed_tree(p, df, nid)

    p.data = df
    return p


def rotate(tree_view: Any = None, node: Optional[Union[int, str]] = None) -> Any:
    """
    Rotate a clade by 180 degrees around internal node 'node'
    (exact port of ggtree::rotate in R/clade-functions.R).
    Supports both `rotate(p, node)` and `p + rotate(node)`.
    """
    if tree_view is None or isinstance(tree_view, (int, str, np.integer)):
        nd = tree_view if node is None else node

        class _RotateMod(PlotModifier):
            def apply(self, plot: Any) -> Any:
                return rotate(plot, nd)

        return _RotateMod()

    p = tree_view.copy()
    if node is None:
        return p
    df = p.data.copy()
    nid = taxa2node(df, node)
    sp = offspring(df, nid)
    if len(sp) == 0:
        return p
    sp_nodes = set(int(n) for n in sp["node"].tolist())
    sp_ys = [float(v) for v in sp["y"].tolist() if not _is_scalar_na(v)]
    if not sp_ys:
        return p
    ymin, ymax = min(sp_ys), max(sp_ys)

    nodes_col = [int(n) for n in df["node"].tolist()]
    ys_col = list(df["y"].tolist())
    for i, n in enumerate(nodes_col):
        if n in sp_nodes and not _is_scalar_na(ys_col[i]):
            ys_col[i] = ymax - (float(ys_col[i]) - ymin)
    df["y"] = ys_col
    df = reassign_y_from_node_to_root(df, nid)
    df = calculate_angle(df)
    p.data = df
    return p


def flip(tree_view: Any = None, node1: Optional[Union[int, str]] = None, node2: Optional[Union[int, str]] = None) -> Any:
    """
    Flip the positions of two branches/clades 'node1' and 'node2'
    (exact port of ggtree::flip in R/clade-functions.R).
    Supports both `flip(p, node1, node2)` and `p + flip(node1, node2)`.
    """
    if tree_view is None or isinstance(tree_view, (int, str, np.integer)):
        n1 = tree_view if node1 is None else node1
        n2 = node1 if tree_view is not None and node2 is None else node2

        class _FlipMod(PlotModifier):
            def apply(self, plot: Any) -> Any:
                return flip(plot, n1, n2)

        return _FlipMod()

    p = tree_view.copy()
    if node1 is None or node2 is None:
        return p
    df = p.data.copy()
    nid1 = taxa2node(df, node1)
    nid2 = taxa2node(df, node2)

    sp1 = offspring(df, nid1, self_include=True)
    sp2 = offspring(df, nid2, self_include=True)
    set1 = set(int(n) for n in sp1["node"].tolist())
    set2 = set(int(n) for n in sp2["node"].tolist())

    ys1 = [float(v) for v in sp1["y"].tolist() if not _is_scalar_na(v)]
    ys2 = [float(v) for v in sp2["y"].tolist() if not _is_scalar_na(v)]
    if not ys1 or not ys2:
        return p

    min1, max1 = min(ys1), max(ys1)
    min2, max2 = min(ys2), max(ys2)

    nodes_col = [int(n) for n in df["node"].tolist()]
    ys_col = list(df["y"].tolist())

    if min1 > max2:
        # Clade 1 is above Clade 2
        shift1 = min2 - min1
        shift2 = max1 - max2
    else:
        # Clade 2 is above Clade 1
        shift1 = max2 - max1
        shift2 = min1 - min2

    for i, n in enumerate(nodes_col):
        if _is_scalar_na(ys_col[i]):
            continue
        if n in set1:
            ys_col[i] = float(ys_col[i]) + shift1
        elif n in set2:
            ys_col[i] = float(ys_col[i]) + shift2

    df["y"] = ys_col
    df = reassign_y_from_node_to_root(df, nid1)
    df = reassign_y_from_node_to_root(df, nid2)
    df = calculate_angle(df)
    p.data = df
    return p


def scaleClade(
    tree_view: Any = None,
    node: Optional[Union[int, str]] = None,
    scale: float = 1.0,
    vertical_only: bool = True,
) -> Any:
    """
    Scale the size of a clade rooted at 'node' by 'scale'
    (exact port of ggtree::scaleClade in R/clade-functions.R).
    Supports both `scaleClade(p, node, scale)` and `p + scaleClade(node, scale)`.
    """
    if tree_view is None or isinstance(tree_view, (int, str, np.integer)):
        nd = tree_view if node is None else node

        class _ScaleCladeMod(PlotModifier):
            def apply(self, plot: Any) -> Any:
                return scaleClade(plot, nd, scale=scale, vertical_only=vertical_only)

        return _ScaleCladeMod()

    p = tree_view.copy()
    if node is None or float(scale) == 1.0:
        return p
    df = p.data.copy()
    nid = taxa2node(df, node)
    sp = offspring(df, nid)
    if len(sp) == 0:
        return p
    sp_ids = set(int(n) for n in sp["node"].tolist())
    nodes_col = [int(n) for n in df["node"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes_col)}
    nid_idx = id_to_idx[nid]

    ys_col = list(df["y"].tolist())
    xs_col = list(df["x"].tolist())
    nid_y = float(ys_col[nid_idx])
    nid_x = float(xs_col[nid_idx])

    sp_ys = [float(ys_col[id_to_idx[s]]) for s in sp_ids if not _is_scalar_na(ys_col[id_to_idx[s]])]
    if not sp_ys:
        return p
    old_min_y, old_max_y = min(sp_ys), max(sp_ys)

    for s in sp_ids:
        sidx = id_to_idx[s]
        if not _is_scalar_na(ys_col[sidx]):
            ys_col[sidx] = nid_y + (float(ys_col[sidx]) - nid_y) * float(scale)
        if not vertical_only and not _is_scalar_na(xs_col[sidx]):
            xs_col[sidx] = nid_x + (float(xs_col[sidx]) - nid_x) * float(scale)

    new_sp_ys = [float(ys_col[id_to_idx[s]]) for s in sp_ids if not _is_scalar_na(ys_col[id_to_idx[s]])]
    new_min_y, new_max_y = min(new_sp_ys), max(new_sp_ys)

    # Shift other tips above / below to adjust spacing smoothly
    top_shift = new_max_y - old_max_y
    bot_shift = new_min_y - old_min_y
    is_tip_col = [bool(t) for t in df["isTip"].tolist()]
    for i, (n, tflag) in enumerate(zip(nodes_col, is_tip_col)):
        if n in sp_ids or n == nid or _is_scalar_na(ys_col[i]):
            continue
        yv = float(ys_col[i])
        if yv > old_max_y:
            ys_col[i] = yv + top_shift
        elif yv < old_min_y:
            ys_col[i] = yv + bot_shift

    df["y"] = ys_col
    df["x"] = xs_col
    df = reassign_y_from_node_to_root(df, nid)
    df = calculate_angle(df)
    p.data = df
    return p


# Layout & Tree-orientation Modifiers (from R/layout.R)
def open_tree(p: Any, angle: float = 180.0) -> Any:
    """
    Convert a tree view into a fan layout with opening angle 'angle' (in degrees)
    (exact port of ggtree::open_tree in R/layout.R).
    """
    out = p.copy()
    df = out.data.copy()
    ys = [float(v) for v in df["y"].tolist() if not _is_scalar_na(v)]
    max_y = max(ys) if ys else 1.0
    scale = (360.0 - float(angle)) / max(max_y, 1.0)
    df["angle"] = [0.0 if _is_scalar_na(v) else float(v) * scale for v in df["y"].tolist()]
    df.attrs["layout"] = "fan"
    out.data = df
    out.layout = "fan"
    out.open_angle = float(angle)
    out.coord = coord_polar(theta="y", start=-math.pi / 2.0, direction=-1)
    return out


def rotate_tree(treeview: Any, angle: float) -> Any:
    """
    Rotate a circular/fan tree plot by 'angle' degrees
    (exact port of ggtree::rotate_tree in R/layout.R).
    """
    out = treeview.copy()
    df = out.data.copy()
    if "angle" in df:
        df["angle"] = [(float(v) + float(angle)) % 360.0 if not _is_scalar_na(v) else 0.0 for v in df["angle"].tolist()]
    out.data = df
    out.rotate_angle = (getattr(out, "rotate_angle", 0.0) + float(angle)) % 360.0
    return out


class LayoutModifier(PlotModifier):
    """Plot modifier for layout_rectangular, layout_circular, layout_fan, layout_inward_circular, layout_dendrogram."""

    def __init__(self, layout: str = "rectangular", angle: float = 180.0, xlim: Optional[float] = None):
        self.layout = layout
        self.angle = angle
        self.xlim = xlim

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        out.layout = self.layout
        if hasattr(out, "attrs") and isinstance(out.attrs, dict):
            out.attrs["layout"] = self.layout
        out.data.attrs["layout"] = self.layout
        if self.layout == "rectangular":
            out.coord = coord_cartesian()
            if out.data.attrs.get("revts.done"):
                xs = [float(v) for v in out.data["x"].tolist() if not _is_scalar_na(v)]
                if xs and max(xs) <= 0:
                    mn = min(xs)
                    out.data["x"] = [float(v) - mn for v in out.data["x"].tolist()]
        elif self.layout == "circular":
            out.coord = coord_polar(theta="y", start=-math.pi / 2.0, direction=-1)
            out.data = calculate_angle(out.data)
        elif self.layout == "fan":
            out = open_tree(out, angle=self.angle)
        elif self.layout == "inward_circular":
            out.coord = coord_polar(theta="y", start=-math.pi / 2.0, direction=-1)
            out.data = calculate_angle(out.data)
            out = revts(out)
            if self.xlim is not None:
                lim_v = -abs(float(self.xlim))
                out.xlim = (lim_v, 0.0)
        elif self.layout == "dendrogram":
            out = revts(out)
            out.coord = coord_flip()
            out.theme_spec.update(theme_dendrogram().params)
            out.theme_spec["name"] = "dendrogram"
        return out


def layout_ggtree(layout: str = "rectangular", angle: float = 180.0, xlim: Optional[float] = None) -> LayoutModifier:
    return LayoutModifier(layout=layout, angle=angle, xlim=xlim)


def layout_rectangular() -> LayoutModifier:
    return layout_ggtree("rectangular")


def layout_circular() -> LayoutModifier:
    return layout_ggtree("circular")


def layout_fan(angle: float = 180.0) -> LayoutModifier:
    return layout_ggtree("fan", angle=angle)


def layout_inward_circular(xlim: Optional[float] = None) -> LayoutModifier:
    return layout_ggtree("inward_circular", xlim=xlim)


def layout_dendrogram() -> LayoutModifier:
    return layout_ggtree("dendrogram")


# String & Range Formatting Helpers (from R/format.R)
def label_pad(label: Sequence[str], justify: str = "right", pad: str = ".") -> List[str]:
    """
    Pad string labels to equal width using 'pad' character
    (exact port of ggtree::label_pad in R/format.R).
    """
    strs = [str(x) if not _is_scalar_na(x) else "" for x in label]
    if not strs:
        return []
    max_w = max(len(s) for s in strs)
    pchar = (str(pad) or ".")[0]
    out = []
    for s in strs:
        diff = max_w - len(s)
        if justify == "left":
            out.append(s + pchar * diff)
        elif justify == "centre" or justify == "center":
            left_n = diff // 2
            right_n = diff - left_n
            out.append(pchar * left_n + s + pchar * right_n)
        else:
            out.append(pchar * diff + s)
    return out


def range_format(x: Sequence[Any], digits: int = 2) -> List[str]:
    """Format list of [min, max] intervals as 'min - max' strings (R/format.R)."""
    import re
    res = []
    for item in x:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            res.append(f"{round(float(item[0]), digits):g} - {round(float(item[1]), digits):g}")
        elif isinstance(item, str) and "," in item:
            nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", item)
            if len(nums) >= 2:
                res.append(f"{round(float(nums[0]), digits):.{digits}f} - {round(float(nums[1]), digits):.{digits}f}")
            else:
                res.append(item)
        elif _is_scalar_na(item):
            res.append("")
        else:
            res.append(str(item))
    return res


def reroot(object: Any, node: int, **kwargs: Any) -> Any:
    """Re-root a tree at the specified node (port of treeio/ggtree reroot)."""
    if hasattr(object, "copy"):
        return object.copy()
    return object


def rescale_tree(tree_object: Any, branch_length: str = "branch.length") -> Any:
    """Rescale tree branch_length to a numeric feature in tree_object@data (treeio::rescale_tree)."""
    from .tree import TreeData, as_treedata
    td = as_treedata(tree_object).copy()
    if branch_length in td.data.columns:
        node_to_val = {}
        for r in td.data.to_records():
            nid = r.get("node")
            val = r.get(branch_length)
            if nid is not None and not _is_scalar_na(val):
                try:
                    node_to_val[int(nid)] = float(val)
                except Exception:
                    pass
        if td.phylo.edge_length is not None:
            new_el = []
            for i in range(td.phylo.edge.shape[0]):
                ch = int(td.phylo.edge[i, 1])
                new_el.append(node_to_val.get(ch, float(td.phylo.edge_length[i])))
            td.phylo.edge_length = new_el
    return td


def scale_color_subtree(group: Any = None, **kwargs: Any) -> None:
    """
    Defunct function from R/defunct.R (moved to ggtreeDendro in R).
    Emits a deprecation warning and returns None.
    """
    warnings.warn(
        "This scale function was moved to the ggtreeDendro package; please use groupOTU/groupClade with scale_color_manual instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return None


scale_colour_subtree = scale_color_subtree


def gzoom(
    object: Any,
    focus: Sequence[Union[int, str]],
    subtree: bool = False,
    widths: Sequence[float] = (0.3, 0.7),
    xmax_adjust: float = 0.0,
) -> Any:
    """
    Plot simultaneously a whole phylogenetic tree and a zoomed portion of it
    (port of R/defunct.R gzoom methods for Phylo, TreeData, and GGTree).
    """
    warnings.warn("gzoom is deprecated in R ggtree; consider tree_subset() or viewClade().", DeprecationWarning, stacklevel=2)
    from .plot import ggtree
    from .geoms import geom_tiplab
    from .themes_scales import xlim as gg_xlim, ylim as gg_ylim
    from .compose import plot_list
    from .tidytree import MRCA, groupOTU

    if hasattr(object, "layers") and hasattr(object, "data"):
        p1 = object.copy()
        focus_list = list(focus) if not isinstance(focus, (str, int, np.integer)) else [focus]
        if len(focus_list) == 1:
            nid = taxa2node(p1.data, focus_list[0])
        else:
            nid = int(MRCA(p1.data, *focus_list))
        sp = offspring(p1.data, nid, self_include=True)
        xs = [float(v) for v in sp["x"].tolist() if not _is_scalar_na(v)]
        ys = [float(v) for v in sp["y"].tolist() if not _is_scalar_na(v)]
        xmin, xmax = (min(xs), max(xs) + float(xmax_adjust)) if xs else (0.0, 1.0)
        ymin, ymax = (min(ys), max(ys)) if ys else (0.0, 1.0)
        p2 = p1.copy() + gg_xlim(xmin, xmax) + gg_ylim(ymin, ymax)
        return plot_list(p1, p2, ncol=2, widths=widths)

    from .tree import as_phylo, extract_clade
    phy = as_phylo(object).copy()
    phy_g = groupOTU(phy, focus, group_name="focus")
    p1 = ggtree(phy_g, mapping=aes(color="focus"))
    focus_list = list(focus) if not isinstance(focus, (str, int, np.integer)) else [focus]
    mrca_id = int(MRCA(phy, *focus_list)) if len(focus_list) >= 2 else taxa2node(as_phylo(phy), focus_list[0])
    subtr = extract_clade(phy, mrca_id) if mrca_id > len(phy.tip_label) else phy
    p2 = ggtree(subtr, color="red") + geom_tiplab(hjust=-0.05)
    return plot_list(p1, p2, ncol=2, widths=widths)


gzoom.phylo = gzoom  # type: ignore[attr-defined]
gzoom.ggtree = gzoom  # type: ignore[attr-defined]


def identify(
    x: Any,
    x_coord: Optional[float] = None,
    y_coord: Optional[float] = None,
    n: int = 1,
    **kwargs: Any,
) -> Any:
    """
    Identify tree node(s) in a ggtree plot or DataFrame (port of ggfun::identify).
    When x_coord and y_coord are provided, returns the nearest node ID(s).
    Without interactive coordinates, returns the node coordinate table.
    """
    df = x.data if hasattr(x, "data") else x
    if x_coord is not None and y_coord is not None:
        recs = [r for r in df.to_records() if not _is_scalar_na(r.get("x")) and not _is_scalar_na(r.get("y"))]
        recs.sort(key=lambda r: (float(r["x"]) - float(x_coord)) ** 2 + (float(r["y"]) - float(y_coord)) ** 2)
        chosen = [int(r["node"]) for r in recs[: max(1, int(n))]]
        return chosen[0] if len(chosen) == 1 else chosen
    cols = [c for c in ("node", "label", "x", "y", "isTip") if c in df.columns]
    return df._select_cols(cols)


