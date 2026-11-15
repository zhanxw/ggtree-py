"""
All ggtree and ggplot2 geometric layers (geoms) and stats:
geom_tree, geom_tree2, geom_tiplab, geom_tiplab2, geom_nodelab, geom_nodelab2,
geom_point2, geom_point, geom_tippoint, geom_nodepoint, geom_rootpoint,
geom_text2, geom_text, geom_label2, geom_label, geom_segment2, geom_segment,
geom_aline, geom_hilight, geom_highlight, geom_hilight_encircle,
geom_cladelab, geom_cladelabel, geom_cladelabel2, geom_strip, geom_striplab,
geom_taxalink, geom_tanglelink, geom_treescale, geom_balance, geom_range,
geom_rootedge, geom_zoom_clade, geom_motif, geom_inset, and standard ggplot2 geoms.
"""

from __future__ import annotations

import math
import re
import warnings
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union
import numpy as np

from .aes import Aes, SymExpr, aes, eval_expr, eval_subset, modifyList
from .dataframe import ColumnSeries, DataFrame, _is_scalar_na, bind_rows, ensure_dataframe, is_dataframe_like
from .layout import (
    adjust_cladelabel_angle,
    check_reverse,
    clade_geometry,
    get_clade_position_,
    get_layout,
)
from .themes_scales import PlotModifier
from .tidytree import child, getRoot, isTip, offspring, taxa2node


def _extract_mapping_and_data(
    mapping: Any = None,
    data: Any = None,
    kwargs: Optional[Dict[str, Any]] = None,
) -> Tuple[Aes, Any]:
    if is_dataframe_like(mapping) or (callable(mapping) and not isinstance(mapping, SymExpr)):
        data, mapping = mapping, data
    if mapping is None:
        m = Aes()
    elif isinstance(mapping, Aes):
        m = mapping.copy()
    elif isinstance(mapping, dict):
        m = Aes(**mapping)
    else:
        m = Aes()
    return m, data


class Layer(PlotModifier):
    """Represents a geometric layer in a ggtree / ggplot object."""

    def __init__(
        self,
        geom_type: str,
        mapping: Optional[Aes] = None,
        data: Any = None,
        params: Optional[Dict[str, Any]] = None,
        stat: str = "identity",
    ):
        self.geom_type = geom_type
        self.mapping = mapping.copy() if isinstance(mapping, Aes) else Aes(**(mapping or {}))
        self.data = data
        self.params: Dict[str, Any] = dict(params) if params else {}
        if "colour" in self.params and "color" not in self.params:
            self.params["color"] = self.params["colour"]
        self.stat = stat

    def copy(self) -> "Layer":
        return Layer(
            geom_type=self.geom_type,
            mapping=self.mapping.copy(),
            data=self.data.copy() if isinstance(self.data, DataFrame) else (ensure_dataframe(self.data) if is_dataframe_like(self.data) else self.data),
            params=dict(self.params),
            stat=self.stat,
        )

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        lyr = self.copy()
        if lyr.params.get("to_bottom", lyr.params.get("to.bottom", False)):
            out.layers.insert(0, lyr)
        else:
            out.layers.append(lyr)
        return out

    def compute_data(self, plot_data: DataFrame, plot_mapping: Aes, layout: str = "rectangular") -> DataFrame:
        """Compute the layer's concrete geometry DataFrame for ggplot_build and rendering."""
        # Resolve layer data
        if self.data is None:
            df = plot_data.copy()
        elif callable(self.data) and not isinstance(self.data, SymExpr):
            df = ensure_dataframe(self.data(plot_data.copy()))
        elif isinstance(self.data, DataFrame):
            df = self.data.copy()
        else:
            df = ensure_dataframe(self.data)

        inherit_flag = self.params.get("inherit.aes", self.params.get("inherit_aes", True))
        merged_aes = plot_mapping.merge(self.mapping) if inherit_flag else self.mapping.copy()
        extra_env = getattr(merged_aes, "_env", None)

        # Apply subset filter if present in mapping or params
        subset_expr = merged_aes.get("subset", self.params.get("subset", None))
        if subset_expr is not None and len(df) > 0:
            mask = eval_subset(df, subset_expr, extra_env=extra_env)
            df = df[mask]

        if self.geom_type in ("tree", "tree_horizontal", "tree_vertical"):
            res = _compute_tree_segments(
                df,
                merged_aes,
                self.params,
                layout=self.params.get("layout") or layout,
                segment_mode=self.geom_type,
            )
            if subset_expr is not None and len(res) > 0:
                res["subset"] = [True] * len(res)
            return res
        if self.geom_type == "rootedge":
            return _compute_rootedge(df, merged_aes, self.params, layout=layout)
        if self.geom_type == "aline":
            res = _compute_aline(df, merged_aes, self.params, plot_data=plot_data)
            if subset_expr is not None and len(res) > 0:
                res["subset"] = [True] * len(res)
            return res
        if self.geom_type == "hilight":
            return _compute_hilight(plot_data, df, merged_aes, self.params, layout=layout)
        if self.geom_type == "balance":
            return _compute_balance(plot_data, self.params, layout=layout)
        if self.geom_type in ("cladelab_bar", "cladelab_text"):
            return _compute_cladelab(plot_data, df, merged_aes, self.params, part=self.geom_type, layout=layout)
        if self.geom_type in ("strip_bar", "strip_text"):
            return _compute_strip(plot_data, df, merged_aes, self.params, part=self.geom_type, layout=layout)
        if self.geom_type == "taxalink":
            return _compute_taxalink(plot_data, df, merged_aes, self.params)
        if self.geom_type == "tanglelink":
            return _compute_tanglelink(plot_data, df, merged_aes, self.params)
        if self.geom_type == "treescale":
            return _compute_treescale(plot_data, self.params)
        if self.geom_type == "range":
            return _compute_range(df, merged_aes, self.params)

        # Standard point / text / label / segment / bar / tile layer evaluation
        out = df.copy()
        for aes_k, aes_expr in merged_aes.items():
            if aes_k == "subset":
                continue
            val = eval_expr(df, aes_expr, extra_env=extra_env)
            if val is not None:
                out[aes_k] = val
                expr_key = str(aes_expr)
                if expr_key in df._levels:
                    out._levels[aes_k] = list(df._levels[expr_key])
                else:
                    out._levels.pop(aes_k, None)

        # Apply offset / nudge adjustments for text/tiplab/nodelab
        offset_val = self.params.get("offset", 0.0)
        nudge_x = self.params.get("nudge_x", 0.0)
        nudge_y = self.params.get("nudge_y", 0.0)
        if self.geom_type in ("tiplab", "nodelab", "text", "label") and "x" in out:
            if self.params.get("align", False) and "x" in plot_data:
                max_x = max((float(v) for v in plot_data["x"].tolist() if not _is_scalar_na(v)), default=0.0)
                out["x"] = [max_x + float(offset_val) + float(nudge_x)] * len(out)
            elif (offset_val or nudge_x) and all(isinstance(v, (int, float, np.number)) for v in out["x"].tolist() if not _is_scalar_na(v)):
                out["x"] = [float(v) + float(offset_val) + float(nudge_x) if not _is_scalar_na(v) else None for v in out["x"].tolist()]
        if nudge_y and "y" in out:
            out["y"] = [float(v) + float(nudge_y) if not _is_scalar_na(v) else None for v in out["y"].tolist()]

        for pk in ("color", "colour", "fill", "size", "alpha", "shape", "linetype", "linewidth", "hjust", "vjust", "angle"):
            if pk in self.params and self.params[pk] is not None and pk not in out:
                out[pk] = [self.params[pk]] * len(out)
        if "color" in out and "colour" not in out:
            out["colour"] = out["color"]
            if "color" in out._levels:
                out._levels["colour"] = list(out._levels["color"])
        elif "colour" in out and "color" not in out:
            out["color"] = out["colour"]
            if "colour" in out._levels:
                out._levels["color"] = list(out._levels["colour"])
        if self.geom_type == "point":
            if "shape" not in out:
                out["shape"] = [19] * len(out)
            if "size" not in out:
                out["size"] = [1.5] * len(out)
            if "colour" not in out:
                out["colour"] = ["black"] * len(out)
                out["color"] = ["black"] * len(out)
        if subset_expr is not None:
            out["subset"] = [True] * len(out)
        if ".panel" in out and "PANEL" not in out:
            out["PANEL"] = out[".panel"]
        elif "PANEL" not in out:
            out["PANEL"] = [self.params.get("panel", "Tree")] * len(out)
        return out


class LayerList(PlotModifier, list):
    """List of layers/modifiers returned by composite geoms (like geom_cladelab, geom_tiplab(align=True))."""

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        for item in self:
            out = out + item
        return out


# Internal geometry builders for specialized ggtree geoms
def _compute_tree_segments(
    df: DataFrame,
    mapping: Aes,
    params: Dict[str, Any],
    layout: str = "rectangular",
    segment_mode: str = "tree",
) -> DataFrame:
    if len(df) == 0:
        return DataFrame(columns=["x", "y", "xend", "yend", "node", "parent"])
    extra_env = getattr(mapping, "_env", None)
    nodes = [int(n) for n in df["node"].tolist()]
    parents = [int(p) for p in df["parent"].tolist()]
    xs: List[Optional[float]] = [None if _is_scalar_na(v) else float(v) for v in df["x"].tolist()]
    ys: List[Optional[float]] = [None if _is_scalar_na(v) else float(v) for v in df["y"].tolist()]

    has_id = ".id" in df
    id_col = [str(v) for v in df[".id"].tolist()] if has_id else None
    if id_col is not None:
        id_to_idx: Dict[Any, int] = {(id_col[i], n): i for i, n in enumerate(nodes)}
    else:
        id_to_idx = {n: i for i, n in enumerate(nodes)}

    color_expr = mapping.get("colour", mapping.get("color"))
    color_levels = df._levels.get(str(color_expr)) if color_expr is not None else None
    if color_expr is not None:
        c_eval = eval_expr(df, color_expr, extra_env=extra_env)
        color_col = c_eval.tolist() if c_eval is not None else [params.get("colour", params.get("color", "black"))] * len(df)
    else:
        color_col = [params.get("colour", params.get("color", "black"))] * len(df)

    def _with_color_levels(res_df: DataFrame) -> DataFrame:
        if color_levels:
            lvl_copy = list(color_levels)
            res_df._levels["color"] = lvl_copy
            res_df._levels["colour"] = list(lvl_copy)
        return res_df

    lw_expr = mapping.get("linewidth", mapping.get("size"))
    default_lw = params.get("linewidth", params.get("size", 0.5))
    if lw_expr is not None:
        lw_eval = eval_expr(df, lw_expr, extra_env=extra_env)
        size_col = lw_eval.tolist() if lw_eval is not None else [default_lw] * len(df)
    else:
        size_col = [default_lw] * len(df)

    if "linetype" in mapping:
        lt_eval = eval_expr(df, mapping["linetype"], extra_env=extra_env)
        ltype_col = lt_eval.tolist() if lt_eval is not None else [params.get("linetype", 1)] * len(df)
    else:
        ltype_col = [params.get("linetype", 1)] * len(df)

    if "alpha" in mapping:
        al_eval = eval_expr(df, mapping["alpha"], extra_env=extra_env)
        alpha_col = al_eval.tolist() if al_eval is not None else [params.get("alpha", 1.0)] * len(df)
    else:
        alpha_col = [params.get("alpha", 1.0)] * len(df)

    direct_layouts = {"slanted", "equal_angle", "daylight", "ape", "tree_and_leaf", "unrooted", "external", "radial"}
    rows: List[Dict[str, Any]] = []

    if segment_mode in ("tree_horizontal", "tree_vertical"):
        if layout in direct_layouts:
            if segment_mode == "tree_vertical":
                return DataFrame(columns=["x", "y", "xend", "yend", "node", "parent"])
            for i, (p, n) in enumerate(zip(parents, nodes)):
                key = (id_col[i], p) if id_col is not None else p
                pidx = id_to_idx.get(key)
                cx, cy = xs[i], ys[i]
                px = xs[pidx] if pidx is not None else cx
                py = ys[pidx] if pidx is not None else cy
                if cx is None or cy is None or px is None or py is None:
                    rx, ry, rxend, ryend = None, None, None, None
                else:
                    rx, ry, rxend, ryend = px, py, cx, cy
                rows.append({
                    "x": rx, "y": ry, "xend": rxend, "yend": ryend,
                    "node": n, "parent": p, "segment_type": "diagonal",
                    "color": color_col[i], "colour": color_col[i],
                    "size": size_col[i], "linewidth": size_col[i],
                    "linetype": ltype_col[i], "alpha": alpha_col[i], "PANEL": "Tree",
                })
            return _with_color_levels(DataFrame(rows))

        # Rectangular-family 2-layer stat behavior (1 row per row in df, matching StatTreeHorizontal/StatTreeVertical)
        for i, (p, n) in enumerate(zip(parents, nodes)):
            key = (id_col[i], p) if id_col is not None else p
            pidx = id_to_idx.get(key)
            cx, cy = xs[i], ys[i]
            px = xs[pidx] if pidx is not None else None
            py = ys[pidx] if pidx is not None else None
            if cx is None or cy is None or px is None or py is None:
                rx, ry, rxend, ryend = None, None, None, None
            elif p == n:
                rx, ry, rxend, ryend = cx, cy, cx, cy
            elif segment_mode == "tree_horizontal":
                rx, ry, rxend, ryend = px, cy, cx, cy
            else:
                rx, ry, rxend, ryend = px, py, px, cy
            seg_kind = "horizontal" if segment_mode == "tree_horizontal" else "vertical"
            rows.append({
                "x": rx, "y": ry, "xend": rxend, "yend": ryend,
                "node": n, "parent": p, "segment_type": seg_kind,
                "color": color_col[i], "colour": color_col[i],
                "size": size_col[i], "linewidth": size_col[i],
                "linetype": ltype_col[i], "alpha": alpha_col[i], "PANEL": "Tree",
            })
        return _with_color_levels(DataFrame(rows))

    for i, (p, n) in enumerate(zip(parents, nodes)):
        key = (id_col[i], p) if id_col is not None else p
        if p == n or key not in id_to_idx:
            continue
        pidx = id_to_idx[key]
        px, py = xs[pidx], ys[pidx]
        cx, cy = xs[i], ys[i]
        if cx is None or cy is None or px is None or py is None:
            continue
        c_val = color_col[i]
        s_val = size_col[i]
        lt_val = ltype_col[i]
        al_val = alpha_col[i]

        if layout in direct_layouts:
            rows.append({
                "x": px, "y": py, "xend": cx, "yend": cy,
                "node": n, "parent": p, "segment_type": "diagonal",
                "color": c_val, "colour": c_val, "size": s_val, "linewidth": s_val,
                "linetype": lt_val, "alpha": al_val, "PANEL": "Tree",
            })
        else:
            # Horizontal edge: (px, cy) -> (cx, cy)
            rows.append({
                "x": px, "y": cy, "xend": cx, "yend": cy,
                "node": n, "parent": p, "segment_type": "horizontal",
                "color": c_val, "colour": c_val, "size": s_val, "linewidth": s_val,
                "linetype": lt_val, "alpha": al_val, "PANEL": "Tree",
            })
            # Vertical connector: (px, py) -> (px, cy)
            rows.append({
                "x": px, "y": py, "xend": px, "yend": cy,
                "node": n, "parent": p, "segment_type": "vertical",
                "color": c_val, "colour": c_val, "size": s_val, "linewidth": s_val,
                "linetype": lt_val, "alpha": al_val, "PANEL": "Tree",
            })
    return _with_color_levels(DataFrame(rows))


def _compute_rootedge(df: DataFrame, mapping: Aes, params: Dict[str, Any], layout: str = "rectangular") -> DataFrame:
    roots = df[df["node"] == df["parent"]]
    if len(roots) == 0:
        return DataFrame(columns=["x", "y", "xend", "yend", "node"])
    r_row = roots.to_records()[0]
    rx = float(r_row.get("x", 0.0) or 0.0)
    ry = float(r_row.get("y", 0.0) or 0.0)
    re_len = params.get("rootedge", None)
    if re_len is None:
        re_len = r_row.get("branch.length", None)
    if _is_scalar_na(re_len):
        return DataFrame(columns=["x", "y", "xend", "yend", "node"])
    d = float(re_len)
    rev = check_reverse(df)
    xend = rx + d if rev else rx - d
    c_val = params.get("color", "black")
    s_val = params.get("linewidth", params.get("size", 0.5))
    return DataFrame([{
        "x": rx, "y": ry, "xend": xend, "yend": ry,
        "node": int(r_row["node"]), "color": c_val, "colour": c_val,
        "size": s_val, "linewidth": s_val, "PANEL": "Tree",
    }])


def _compute_aline(
    df: DataFrame,
    mapping: Aes,
    params: Dict[str, Any],
    plot_data: Optional[DataFrame] = None,
) -> DataFrame:
    tips = df[df["isTip"]] if "isTip" in df else df
    if len(tips) == 0:
        return DataFrame(columns=["x", "y", "xend", "yend", "label"])
    ref_df = plot_data if (plot_data is not None and "x" in plot_data) else tips
    if "isTip" in ref_df:
        ref_tips = ref_df[ref_df["isTip"]]
        if len(ref_tips) > 0:
            ref_df = ref_tips
    xs = [float(v) for v in ref_df["x"].tolist() if not _is_scalar_na(v)]
    max_x = (max(xs) if xs else 0.0) + float(params.get("offset", 0.0))
    col_val = params.get("colour", params.get("color", "black"))
    lw_val = params.get("linewidth", params.get("size", 0.5))
    rows = []
    for r in tips.to_records():
        if _is_scalar_na(r.get("x")) or _is_scalar_na(r.get("y")):
            continue
        rx = float(r["x"])
        ry = float(r["y"])
        rows.append({
            "x": rx,
            "xend": max_x,
            "y": ry,
            "yend": ry,
            "node": r.get("node"),
            "label": r.get("label"),
            "linetype": params.get("linetype", "dotted"),
            "size": lw_val,
            "linewidth": lw_val,
            "color": col_val,
            "colour": col_val,
            "PANEL": "Tree",
        })
    return DataFrame(rows)


def _compute_hilight(
    plot_data: DataFrame,
    layer_df: DataFrame,
    mapping: Aes,
    params: Dict[str, Any],
    layout: str = "rectangular",
) -> DataFrame:
    node_param = params.get("node", None)
    extend = params.get("extend", 0.0)
    extendto = params.get("extendto", None)
    align = str(params.get("align", "none"))
    fill_default = params.get("fill", "steelblue")
    alpha_default = params.get("alpha", 0.5)
    color_default = params.get("color", None)
    rev = check_reverse(plot_data)

    extra_env = getattr(mapping, "_env", None)
    if node_param is not None and "node" not in mapping:
        node_list = [node_param] if isinstance(node_param, (int, str, np.integer)) else list(node_param)
        fill_list = [fill_default] * len(node_list)
        alpha_list = [alpha_default] * len(node_list)
    else:
        node_col = eval_expr(layer_df, mapping.get("node", "node"), extra_env=extra_env)
        node_list = node_col.tolist() if node_col is not None else []
        fill_col = eval_expr(layer_df, mapping["fill"], extra_env=extra_env) if "fill" in mapping else None
        fill_list = fill_col.tolist() if fill_col is not None else [fill_default] * len(node_list)
        alpha_col = eval_expr(layer_df, mapping["alpha"], extra_env=extra_env) if "alpha" in mapping else None
        alpha_list = alpha_col.tolist() if alpha_col is not None else [alpha_default] * len(node_list)

    all_xs = [float(v) for v in plot_data["x"].tolist() if not _is_scalar_na(v)]
    tree_min_x = min(all_xs) if all_xs else 0.0
    tree_max_x = max(all_xs) if all_xs else 1.0

    rows: List[Dict[str, Any]] = []
    for nid_raw, f_val, a_val in zip(node_list, fill_list, alpha_list):
        if _is_scalar_na(nid_raw):
            continue
        pos = get_clade_position_(plot_data, nid_raw, reverse=rev).to_records()[0]
        xmin, xmax, ymin, ymax = pos["xmin"], pos["xmax"], pos["ymin"], pos["ymax"]

        if isinstance(extend, (list, tuple)) and len(extend) >= 2:
            xmin -= float(extend[0])
            xmax += float(extend[1])
        else:
            xmax += float(extend)

        if extendto is not None and not _is_scalar_na(extendto):
            xmax = float(extendto)

        if align == "left":
            xmin = tree_min_x
        elif align == "right":
            xmax = tree_max_x
        elif align == "both":
            xmin = tree_min_x
            xmax = tree_max_x

        rows.append({
            "xmin": xmin,
            "xmax": xmax,
            "ymin": ymin,
            "ymax": ymax,
            "clade_root_node": taxa2node(plot_data, nid_raw),
            "fill": f_val,
            "alpha": a_val,
            "color": color_default,
            "colour": color_default,
            "type": params.get("type", "rect"),
            "PANEL": "Tree",
        })
    return DataFrame(rows)


def _compute_balance(plot_data: DataFrame, params: Dict[str, Any], layout: str = "rectangular") -> DataFrame:
    node_val = params.get("node")
    nid = taxa2node(plot_data, node_val)
    if isTip(plot_data, nid):
        raise ValueError(f"Selected node {nid} is a tip. Please select an internal node for geom_balance.")
    ch_df = child(plot_data, nid)
    ch_nodes = [int(n) for n in ch_df["node"].tolist()]
    if len(ch_nodes) < 2:
        warnings.warn(f"Balance node {nid} has fewer than 2 descendants; skipping geom_balance.")
        return DataFrame(columns=["xmin", "xmax", "ymin", "ymax", "node"])
    if len(ch_nodes) > 2:
        warnings.warn(f"Balance node {nid} has more than 2 direct children (polytomy); skipping geom_balance.")
        return DataFrame(columns=["xmin", "xmax", "ymin", "ymax", "node"])

    rev = check_reverse(plot_data)
    extend = float(params.get("extend", 0.0))
    extendto = params.get("extendto", None)
    fill_val = params.get("fill", "steelblue")
    color_val = params.get("color", "white")
    alpha_val = params.get("alpha", 0.5)

    rows = []
    for idx, cid in enumerate(ch_nodes):
        pos = get_clade_position_(plot_data, cid, reverse=rev).to_records()[0]
        xmax = pos["xmax"] + extend
        if extendto is not None and not _is_scalar_na(extendto):
            xmax = float(extendto)
        c_fill = fill_val[idx % len(fill_val)] if isinstance(fill_val, (list, tuple)) and len(fill_val) > 0 else fill_val
        c_color = color_val[idx % len(color_val)] if isinstance(color_val, (list, tuple)) and len(color_val) > 0 else color_val
        rows.append({
            "xmin": pos["xmin"],
            "xmax": xmax,
            "ymin": pos["ymin"],
            "ymax": pos["ymax"],
            "node": cid,
            "parent_balance_node": nid,
            "fill": c_fill,
            "color": c_color,
            "colour": c_color,
            "alpha": alpha_val,
            "PANEL": "Tree",
        })
    return DataFrame(rows)


def _compute_cladelab(
    plot_data: DataFrame,
    layer_df: DataFrame,
    mapping: Aes,
    params: Dict[str, Any],
    part: str = "cladelab_bar",
    layout: str = "rectangular",
) -> DataFrame:
    node_param = params.get("node", None)
    label_param = params.get("label", "")
    extra_env = getattr(mapping, "_env", None)

    raw_col = params.get("colour", params.get("color", None))
    if isinstance(raw_col, (list, tuple, np.ndarray)) and not isinstance(raw_col, (str, bytes)):
        col_seq = list(raw_col)
        def_bcol = col_seq[0] if len(col_seq) >= 1 else "black"
        def_tcol = col_seq[1] if len(col_seq) >= 2 else def_bcol
    else:
        def_bcol = raw_col or "black"
        def_tcol = raw_col or "black"
    default_col = params.get("barcolor") or def_bcol if part == "cladelab_bar" else params.get("textcolor") or def_tcol

    if node_param is not None and "node" not in mapping:
        node_list = [node_param] if isinstance(node_param, (int, str, np.integer)) else list(node_param)
        label_list = [label_param] * len(node_list) if isinstance(label_param, (str, int, float)) else list(label_param)
        offset_list = [float(params.get("offset", 0.0))] * len(node_list)
        color_list = [default_col] * len(node_list)
    else:
        node_col = eval_expr(layer_df, mapping.get("node", "node"), extra_env=extra_env)
        node_list = node_col.tolist() if node_col is not None else []
        lbl_col = eval_expr(layer_df, mapping.get("label", "label"), extra_env=extra_env)
        label_list = lbl_col.tolist() if lbl_col is not None else [label_param] * len(node_list)
        off_col = eval_expr(layer_df, mapping["offset"], extra_env=extra_env) if "offset" in mapping else None
        offset_list = [float(v) for v in off_col.tolist()] if off_col is not None else [float(params.get("offset", 0.0))] * len(node_list)
        c_expr = mapping.get("colour", mapping.get("color"))
        col_col = eval_expr(layer_df, c_expr, extra_env=extra_env) if c_expr is not None else None
        color_list = col_col.tolist() if col_col is not None else [default_col] * len(node_list)

    if not node_list:
        return DataFrame()

    offset_text = float(params.get("offset.text", params.get("offset_text", 0.0)))
    total_offsets = [o + (offset_text if part == "cladelab_text" else 0.0) for o in offset_list]

    geom_res = clade_geometry(
        trdf=plot_data,
        nodeids=node_list,
        offset=total_offsets,
        align=params.get("align", False),
        extend=params.get("extend", 0.0),
        angle=params.get("angle", 0.0),
        horizontal=params.get("horizontal", True),
        label=label_list,
    )
    out = geom_res["bar"] if part == "cladelab_bar" else geom_res["label"]
    out["color"] = color_list
    out["colour"] = color_list
    out["parse"] = [bool(params.get("parse", False))] * len(out)
    if part == "cladelab_text":
        out["size"] = [float(params.get("fontsize", 3.88))] * len(out)
        out["hjust"] = [float(params.get("hjust", 0.0))] * len(out)
        out["family"] = [str(params.get("family", "sans"))] * len(out)
    else:
        bs = float(params.get("barsize", 0.5))
        out["size"] = [bs] * len(out)
        out["linewidth"] = [bs] * len(out)
    out["PANEL"] = ["Tree"] * len(out)
    return out


def _compute_strip(
    plot_data: DataFrame,
    layer_df: DataFrame,
    mapping: Aes,
    params: Dict[str, Any],
    part: str = "strip_bar",
    layout: str = "rectangular",
) -> DataFrame:
    t1_param = params.get("taxa1", None)
    t2_param = params.get("taxa2", None)
    lbl_param = params.get("label", None)
    extra_env = getattr(mapping, "_env", None)

    if t1_param is not None and t2_param is not None and "taxa1" not in mapping:
        t1_list = [t1_param] if isinstance(t1_param, (int, str, np.integer)) else list(t1_param)
        t2_list = [t2_param] if isinstance(t2_param, (int, str, np.integer)) else list(t2_param)
        lbl_list = [lbl_param] * len(t1_list) if isinstance(lbl_param, (str, int, float)) or lbl_param is None else list(lbl_param)
        color_list = [params.get("colour", params.get("color", "black"))] * len(t1_list)
    else:
        t1_col = eval_expr(layer_df, mapping.get("taxa1", "taxa1"), extra_env=extra_env)
        t2_col = eval_expr(layer_df, mapping.get("taxa2", "taxa2"), extra_env=extra_env)
        t1_list = t1_col.tolist() if t1_col is not None else []
        t2_list = t2_col.tolist() if t2_col is not None else []
        lbl_col = eval_expr(layer_df, mapping["label"], extra_env=extra_env) if "label" in mapping else None
        lbl_list = lbl_col.tolist() if lbl_col is not None else [lbl_param] * len(t1_list)
        c_expr = mapping.get("colour", mapping.get("color"))
        col_col = eval_expr(layer_df, c_expr, extra_env=extra_env) if c_expr is not None else None
        color_list = col_col.tolist() if col_col is not None else [params.get("colour", params.get("color", "black"))] * len(t1_list)

    align = bool(params.get("align", True))
    offset = float(params.get("offset", 0.0))
    offset_text = float(params.get("offset.text", params.get("offset_text", 0.0)))
    extend_raw = params.get("extend", 0.0)
    ext = (float(extend_raw[0]), float(extend_raw[1])) if isinstance(extend_raw, (list, tuple)) and len(extend_raw) >= 2 else (float(extend_raw), float(extend_raw))
    raw_angle = params.get("angle", 0.0)
    horiz = bool(params.get("horizontal", True))
    parse_val = bool(params.get("parse", False))

    max_tree_x = max((float(v) for v in plot_data["x"].tolist() if not _is_scalar_na(v)), default=0.0)
    rows = []
    for t1, t2, lbl, col in zip(t1_list, t2_list, lbl_list, color_list):
        n1 = taxa2node(plot_data, t1)
        n2 = taxa2node(plot_data, t2)
        r1 = plot_data[plot_data["node"] == n1].to_records()[0]
        r2 = plot_data[plot_data["node"] == n2].to_records()[0]
        xx = [float(r1["x"]), float(r2["x"])]
        yy = [float(r1["y"]), float(r2["y"])]
        mx = max_tree_x if align else max(xx)
        y_min = min(yy) - ext[1]
        y_max = max(yy) + ext[0]
        ratio = 1.02 if part == "strip_bar" else 1.03
        extra = offset + (offset_text if part == "strip_text" else 0.0)
        x_pos = mx * ratio + extra
        if raw_angle == "auto":
            if layout in ("circular", "fan", "inward_circular", "radial"):
                a1 = float(r1.get("angle", 0.0) or 0.0)
                a2 = float(r2.get("angle", 0.0) or 0.0)
                ang_val = adjust_cladelabel_angle((a1 + a2) / 2.0, horizontal=horiz)
            else:
                ang_val = 0.0
        else:
            ang_val = float(raw_angle)
        rows.append({
            "x": x_pos,
            "xend": x_pos,
            "y": y_min if part == "strip_bar" else (y_min + y_max) / 2.0,
            "yend": y_max,
            "taxa1": t1,
            "taxa2": t2,
            "label": lbl if lbl is not None else "",
            "color": col,
            "colour": col,
            "angle": ang_val,
            "parse": parse_val,
            "fontsize": float(params.get("fontsize", 3.88)),
            "size": float(params.get("barsize", 0.5)) if part == "strip_bar" else float(params.get("fontsize", 3.88)),
            "hjust": float(params.get("hjust", 0.0)),
            "PANEL": "Tree",
        })
    return DataFrame(rows)


def _compute_taxalink(plot_data: DataFrame, layer_df: DataFrame, mapping: Aes, params: Dict[str, Any]) -> DataFrame:
    t1_param = params.get("taxa1", None)
    t2_param = params.get("taxa2", None)
    if t1_param is not None and t2_param is not None and "taxa1" not in mapping:
        t1_list = [t1_param] if isinstance(t1_param, (int, str, np.integer)) else list(t1_param)
        t2_list = [t2_param] if isinstance(t2_param, (int, str, np.integer)) else list(t2_param)
        color_list = [params.get("color", "black")] * len(t1_list)
    else:
        t1_col = eval_expr(layer_df, mapping.get("taxa1", layer_df.columns[0] if len(layer_df.columns) >= 1 else "taxa1"))
        t2_col = eval_expr(layer_df, mapping.get("taxa2", layer_df.columns[1] if len(layer_df.columns) >= 2 else "taxa2"))
        t1_list = t1_col.tolist() if t1_col is not None else []
        t2_list = t2_col.tolist() if t2_col is not None else []
        col_col = eval_expr(layer_df, mapping["color"]) if "color" in mapping else None
        color_list = col_col.tolist() if col_col is not None else [params.get("color", "black")] * len(t1_list)

    offset = float(params.get("offset", 0.0) or 0.0)
    rows = []
    for t1, t2, col in zip(t1_list, t2_list, color_list):
        n1 = taxa2node(plot_data, t1)
        n2 = taxa2node(plot_data, t2)
        r1 = plot_data[plot_data["node"] == n1].to_records()[0]
        r2 = plot_data[plot_data["node"] == n2].to_records()[0]
        rows.append({
            "x": float(r1["x"]) + offset,
            "y": float(r1["y"]),
            "xend": float(r2["x"]) + offset,
            "yend": float(r2["y"]),
            "node1": n1,
            "node2": n2,
            "taxa1": t1,
            "taxa2": t2,
            "curvature": float(params.get("curvature", 0.5)),
            "color": col,
            "colour": col,
            "alpha": float(params.get("alpha", 1.0)),
            "linetype": params.get("linetype", 1),
            "size": float(params.get("size", 0.5)),
            "PANEL": "Tree",
        })
    return DataFrame(rows)


def _compute_tanglelink(plot_data: DataFrame, layer_df: DataFrame, mapping: Aes, params: Dict[str, Any]) -> DataFrame:
    """Exact port of R/method-ggplot-add-tanglegram.R for geom_tanglelink."""
    extra_env = getattr(mapping, "_env", None)
    assoc = params.get("assoc", None)
    if assoc is None:
        assoc = plot_data.attrs.get("tangle_assoc", None)
    if assoc is None:
        raise ValueError("geom_tanglelink requires an 'assoc' data.frame (either passed directly or stored in the plot via ggdoubletree / as_tanglegram).")
    if not isinstance(assoc, DataFrame):
        assoc = DataFrame(assoc)
    if len(assoc.columns) < 2:
        raise ValueError("'assoc' must have at least 2 columns mapping left tips to right tips.")

    curv_val = float(params.get("curvature", 0.0))
    alp_val = float(params.get("alpha", 0.8))
    lt_val = params.get("linetype", 1)
    lw_val = float(params.get("linewidth", params.get("size", 0.5)))

    if all(c in assoc.columns for c in ("x", "y", "xend", "yend")):
        link_df = assoc.copy()
        n = len(link_df)
        if "curvature" not in link_df:
            link_df["curvature"] = [curv_val] * n
        if "alpha" not in link_df:
            link_df["alpha"] = [alp_val] * n
        if "linetype" not in link_df:
            link_df["linetype"] = [lt_val] * n
        if "linewidth" not in link_df:
            link_df["linewidth"] = [lw_val] * n
        if "PANEL" not in link_df:
            link_df["PANEL"] = ["Tree"] * n
    else:
        tips = plot_data[plot_data["isTip"]] if "isTip" in plot_data else plot_data
        if "side" in tips:
            left_tips = tips[tips["side"] == "left"]
            right_tips = tips[tips["side"] == "right"]
        elif "tree" in tips:
            tree_vals = tips["tree"].tolist()
            uniq_trees = []
            for tv in tree_vals:
                if tv not in uniq_trees:
                    uniq_trees.append(tv)
            left_tag = uniq_trees[0] if len(uniq_trees) >= 1 else "left"
            right_tag = uniq_trees[1] if len(uniq_trees) >= 2 else "right"
            left_tips = tips[tips["tree"] == left_tag]
            right_tips = tips[tips["tree"] == right_tag]
        else:
            left_tips = tips
            right_tips = tips

        left_map = {str(r["label"]): (float(r["x"]), float(r["y"])) for r in left_tips.to_records() if not _is_scalar_na(r.get("label"))}
        right_map = {str(r["label"]): (float(r["x"]), float(r["y"])) for r in right_tips.to_records() if not _is_scalar_na(r.get("label"))}

        c0, c1 = assoc.columns[0], assoc.columns[1]
        extra_cols = assoc.columns[2:]

        link_rows: List[Dict[str, Any]] = []
        for r in assoc.to_records():
            l_lbl = str(r.get(c0, ""))
            r_lbl = str(r.get(c1, ""))
            if l_lbl in left_map and r_lbl in right_map:
                lx, ly = left_map[l_lbl]
                rx, ry = right_map[r_lbl]
                rec: Dict[str, Any] = {
                    "x": lx,
                    "y": ly,
                    "xend": rx,
                    "yend": ry,
                    "left": l_lbl,
                    "right": r_lbl,
                    "curvature": curv_val,
                    "alpha": alp_val,
                    "linetype": lt_val,
                    "linewidth": lw_val,
                    "PANEL": "Tree",
                }
                for ec in extra_cols:
                    rec[ec] = r.get(ec)
                link_rows.append(rec)

        if not link_rows:
            warnings.warn("No matching tip pairs found between left and right trees in geom_tanglelink.")
            return DataFrame(columns=["x", "y", "xend", "yend", "left", "right", "PANEL"])
        link_df = DataFrame(link_rows)

    if "subset" in mapping and len(link_df) > 0:
        mask = eval_subset(link_df, mapping["subset"], extra_env=extra_env)
        link_df = link_df[mask]

    color_expr = mapping.get("colour", mapping.get("color"))
    if color_expr is not None:
        c_col = eval_expr(link_df, color_expr, extra_env=extra_env)
        link_df["color"] = c_col
        link_df["colour"] = c_col
    else:
        c_def = params.get("colour", params.get("color", "grey50"))
        link_df["color"] = [c_def] * len(link_df)
        link_df["colour"] = [c_def] * len(link_df)

    for aes_k in ("alpha", "linetype", "linewidth", "size"):
        if aes_k in mapping:
            v_col = eval_expr(link_df, mapping[aes_k], extra_env=extra_env)
            if v_col is not None:
                link_df[aes_k] = v_col
    return link_df


def _compute_treescale(plot_data: DataFrame, params: Dict[str, Any]) -> DataFrame:
    xs = [float(v) for v in plot_data["x"].tolist() if not _is_scalar_na(v)]
    ys = [float(v) for v in plot_data["y"].tolist() if not _is_scalar_na(v)]
    min_x, max_x = (min(xs), max(xs)) if xs else (0.0, 1.0)
    min_y, max_y = (min(ys), max(ys)) if ys else (0.0, 1.0)
    dx = max(max_x - min_x, 1e-6)
    dy = max(max_y - min_y, 1.0)

    x_pos = float(params["x"]) if params.get("x") is not None else min_x + dx * 0.05
    y_pos = float(params["y"]) if params.get("y") is not None else min_y
    width = float(params["width"]) if params.get("width") is not None else round(dx * 0.15, 2) or 0.1
    offset = float(params["offset"]) if params.get("offset") is not None else dy * 0.04

    return DataFrame([{
        "x": x_pos,
        "xend": x_pos + width,
        "y": y_pos,
        "yend": y_pos,
        "label_x": x_pos + width / 2.0,
        "label_y": y_pos + offset,
        "label": f"{width:g}",
        "color": params.get("color", "black"),
        "colour": params.get("color", "black"),
        "size": params.get("linesize", 0.5),
        "fontsize": params.get("fontsize", 3.88),
        "PANEL": "Tree",
    }])


def _compute_range(df: DataFrame, mapping: Aes, params: Dict[str, Any]) -> DataFrame:
    range_col_name = params.get("range", "length_0.95_HPD")
    center_col_name = params.get("center", "auto")
    if range_col_name not in df:
        m_pct = re.search(r"_0\.(\d+)_", range_col_name)
        if m_pct:
            cand = range_col_name.replace(f"_0.{m_pct.group(1)}_", f"_{m_pct.group(1)}%_")
            if cand in df:
                range_col_name = cand
        elif "%" in range_col_name:
            m_p = re.search(r"_(\d+)%_", range_col_name)
            if m_p:
                cand = range_col_name.replace(f"_{m_p.group(1)}%_", f"_0.{m_p.group(1)}_")
                if cand in df:
                    range_col_name = cand
    if range_col_name not in df:
        return DataFrame(columns=["x", "xend", "y", "yend", "node"])
    rows = []
    for r in df.to_records():
        rng = r.get(range_col_name)
        if not isinstance(rng, (list, tuple)) or len(rng) < 2:
            continue
        lo, hi = float(rng[0]), float(rng[1])
        if not (math.isfinite(lo) and math.isfinite(hi)) or abs(hi - lo) <= 1e-9:
            continue
        lo_v, hi_v = min(lo, hi), max(lo, hi)
        node_x = float(r.get("x", 0.0) or 0.0)
        if center_col_name in ("x", "branch"):
            cx = float(r.get(center_col_name, node_x) or node_x)
            half = (hi_v - lo_v) / 2.0
            x_lo, x_hi = cx - half, cx + half
        elif center_col_name != "auto" and center_col_name in df and not _is_scalar_na(r.get(center_col_name)):
            cval = float(r[center_col_name])
            x_lo = cval + node_x - hi_v
            x_hi = cval + node_x - lo_v
        else:
            cval = (lo_v + hi_v) / 2.0
            x_lo = cval + node_x - hi_v
            x_hi = cval + node_x - lo_v
        ry = float(r.get("y", 0.0) or 0.0)
        rows.append({
            "x": x_lo,
            "xend": x_hi,
            "y": ry,
            "yend": ry,
            "node": r.get("node"),
            "color": params.get("color", "red"),
            "colour": params.get("color", "red"),
            "size": params.get("size", 2.0),
            "alpha": params.get("alpha", 0.6),
            "PANEL": "Tree",
        })
    return DataFrame(rows)


# Public Geom Functions (R-compatible API)
def geom_tree(
    mapping: Any = None,
    data: Any = None,
    layout: Optional[str] = None,
    multiPhylo: bool = False,
    continuous: str = "none",
    position: str = "identity",
    **kwargs: Any,
) -> Union[Layer, LayerList]:
    """Add tree branches layer (R/geom_tree.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    params = {"layout": layout, "multiPhylo": multiPhylo, "continuous": continuous, "position": position, **kwargs}
    direct_layouts = {"slanted", "equal_angle", "daylight", "ape", "tree_and_leaf", "unrooted", "external", "radial"}
    if layout is not None and str(layout) in direct_layouts:
        return Layer("tree", mapping=m, data=d, params=params, stat="tree")
    return LayerList([
        Layer("tree_horizontal", mapping=m.copy(), data=d, params=dict(params), stat="tree_horizontal"),
        Layer("tree_vertical", mapping=m.copy(), data=d, params=dict(params), stat="tree_vertical"),
    ])


def geom_tree2(mapping: Any = None, data: Any = None, **kwargs: Any) -> Union[Layer, LayerList]:
    return geom_tree(mapping=mapping, data=data, **kwargs)


def stat_tree(mapping: Any = None, data: Any = None, **kwargs: Any) -> Union[Layer, LayerList]:
    return geom_tree(mapping=mapping, data=data, **kwargs)


def geom_tiplab(
    mapping: Any = None,
    hjust: float = 0.0,
    offset: float = 0.0,
    align: bool = False,
    linetype: Any = "dotted",
    linesize: float = 0.5,
    geom: str = "text",
    as_ylab: bool = False,
    data: Any = None,
    **kwargs: Any,
) -> Union[Layer, LayerList]:
    """Add tip labels layer (R/geom_tiplab.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    saved_env = dict(getattr(m, "_env", {}))
    if "label" not in m:
        m["label"] = "label"
    if "subset" in m:
        orig_sub = m["subset"]
        m["subset"] = SymExpr(
            f"(isTip & ({orig_sub}))",
            lambda df, _s=orig_sub, _e=saved_env: (df["isTip"] == True) & eval_subset(df, _s, extra_env=_e),
        )
    else:
        m["subset"] = SymExpr("isTip", lambda df: df["isTip"] == True)
    m._env = saved_env

    params = {
        "hjust": hjust,
        "offset": offset,
        "align": align,
        "linetype": linetype,
        "linesize": linesize,
        "geom": geom,
        "as_ylab": as_ylab,
        **kwargs,
    }
    text_layer = Layer("tiplab", mapping=m, data=d, params=params)
    if align:
        aline_mapping = Aes(subset=m["subset"], _env=saved_env)
        aline_layer = geom_aline(
            mapping=aline_mapping,
            linetype=linetype,
            size=linesize,
            offset=offset,
            data=d,
        )
        return LayerList([aline_layer, text_layer])
    return text_layer


def geom_tiplab2(mapping: Any = None, hjust: float = 0.0, **kwargs: Any) -> Union[Layer, LayerList]:
    """Add tip labels for circular layouts (R/geom_tiplab.R)."""
    return geom_tiplab(mapping=mapping, hjust=hjust, circular_auto=True, **kwargs)


def geom_nodelab(
    mapping: Any = None,
    nudge_x: float = 0.0,
    nudge_y: float = 0.0,
    geom: str = "text",
    hjust: float = 0.5,
    node: str = "internal",
    data: Any = None,
    **kwargs: Any,
) -> Layer:
    """Add internal node labels layer (R/geom_nodelab.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    saved_env = dict(getattr(m, "_env", {}))
    if "label" not in m:
        m["label"] = "label"
    if "subset" in m:
        orig_sub = m["subset"]
        if node == "internal":
            m["subset"] = SymExpr(
                f"(!isTip & ({orig_sub}))",
                lambda df, _s=orig_sub, _e=saved_env: (df["isTip"] == False) & eval_subset(df, _s, extra_env=_e),
            )
        elif node == "external":
            m["subset"] = SymExpr(
                f"(isTip & ({orig_sub}))",
                lambda df, _s=orig_sub, _e=saved_env: (df["isTip"] == True) & eval_subset(df, _s, extra_env=_e),
            )
    else:
        if node == "internal":
            m["subset"] = SymExpr("!isTip", lambda df: df["isTip"] == False)
        elif node == "external":
            m["subset"] = SymExpr("isTip", lambda df: df["isTip"] == True)
    m._env = saved_env
    params = {"nudge_x": nudge_x, "nudge_y": nudge_y, "geom": geom, "hjust": hjust, "node": node, **kwargs}
    return Layer("nodelab", mapping=m, data=d, params=params)


def geom_nodelab2(mapping: Any = None, **kwargs: Any) -> Layer:
    return geom_nodelab(mapping=mapping, circular_auto=True, **kwargs)


def geom_point2(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    """ggplot2::geom_point with subset support (R/geom_point.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("point", mapping=m, data=d, params=kwargs)


def geom_point(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    return geom_point2(mapping=mapping, data=data, **kwargs)


def geom_tippoint(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    """Add points at tree tips (R/geom_point.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    saved_env = dict(getattr(m, "_env", {}))
    if "subset" in m:
        orig_sub = m["subset"]
        m["subset"] = SymExpr(
            f"(isTip & ({orig_sub}))",
            lambda df, _s=orig_sub, _e=saved_env: (df["isTip"] == True) & eval_subset(df, _s, extra_env=_e),
        )
    else:
        m["subset"] = SymExpr("isTip", lambda df: df["isTip"] == True)
    m._env = saved_env
    return Layer("point", mapping=m, data=d, params=kwargs)


def geom_nodepoint(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    """Add points at internal nodes (R/geom_point.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    saved_env = dict(getattr(m, "_env", {}))
    if "subset" in m:
        orig_sub = m["subset"]
        m["subset"] = SymExpr(
            f"(!isTip & ({orig_sub}))",
            lambda df, _s=orig_sub, _e=saved_env: (df["isTip"] == False) & eval_subset(df, _s, extra_env=_e),
        )
    else:
        m["subset"] = SymExpr("!isTip", lambda df: df["isTip"] == False)
    m._env = saved_env
    return Layer("point", mapping=m, data=d, params=kwargs)


def geom_rootpoint(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    """Add point at root node (R/geom_point.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    saved_env = dict(getattr(m, "_env", {}))
    m["subset"] = SymExpr("node == parent", lambda df: df["node"] == df["parent"])
    m._env = saved_env
    return Layer("point", mapping=m, data=d, params=kwargs)


def geom_text2(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    """ggplot2::geom_text with subset support (R/geom_text.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("text", mapping=m, data=d, params=kwargs)


def geom_text(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    return geom_text2(mapping=mapping, data=data, **kwargs)


def geom_label2(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    """ggplot2::geom_label with subset support (R/geom_label.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("label", mapping=m, data=d, params=kwargs)


def geom_label(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    return geom_label2(mapping=mapping, data=data, **kwargs)


def geom_segment2(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    """ggplot2::geom_segment with subset support (R/geom_segment.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("segment", mapping=m, data=d, params=kwargs)


def geom_segment(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    return geom_segment2(mapping=mapping, data=data, **kwargs)


def geom_aline(mapping: Any = None, linetype: Any = "dotted", size: float = 0.5, data: Any = None, **kwargs: Any) -> Layer:
    """Add horizontal leader lines from tips to aligned tip labels (R/geom_segment.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("aline", mapping=m, data=d, params={"linetype": linetype, "size": size, **kwargs})


def geom_hilight(
    data: Any = None,
    mapping: Any = None,
    node: Any = None,
    fill: str = "steelblue",
    alpha: float = 0.5,
    extend: Any = 0.0,
    extendto: Optional[float] = None,
    align: str = "none",
    type: str = "auto",
    gradient_direction: str = "rt",
    to_bottom: bool = False,
    **kwargs: Any,
) -> Layer:
    """Highlight clade(s) with rectangle, roundrect, gradient, or encircle (R/geom_hilight.R)."""
    if isinstance(data, (int, np.integer)) and node is None:
        node, data = data, None
    m, d = _extract_mapping_and_data(mapping, data)
    gradient_direction = kwargs.pop("gradient.direction", gradient_direction)
    to_bottom = kwargs.pop("to.bottom", to_bottom)
    params = {
        "node": node,
        "fill": fill,
        "alpha": alpha,
        "extend": extend,
        "extendto": extendto,
        "align": align,
        "type": type,
        "gradient_direction": gradient_direction,
        "to_bottom": to_bottom,
        **kwargs,
    }
    return Layer("hilight", mapping=m, data=d, params=params)


geom_highlight = geom_hilight


def geom_hilight_encircle(data: Any = None, mapping: Any = None, node: Any = None, **kwargs: Any) -> Layer:
    return geom_hilight(data=data, mapping=mapping, node=node, type="encircle", **kwargs)


def geom_cladelab(
    node: Any = None,
    label: Any = None,
    data: Any = None,
    mapping: Any = None,
    offset: float = 0.0,
    offset_text: float = 0.0,
    extend: Any = 0.0,
    align: bool = False,
    barsize: float = 0.5,
    barcolour: Optional[str] = None,
    barcolor: Optional[str] = None,
    textcolor: Optional[str] = None,
    textcolour: Optional[str] = None,
    angle: Any = 0.0,
    fontsize: float = 3.88,
    horizontal: bool = True,
    geom: str = "text",
    hjust: float = 0.0,
    **kwargs: Any,
) -> LayerList:
    """Annotate clade(s) with a vertical bar and text/label/image (R/geom_cladelab.R)."""
    offset_text = kwargs.pop("offset.text", offset_text)
    m, d = _extract_mapping_and_data(mapping, data)
    bcol = barcolor or barcolour or kwargs.get("color", kwargs.get("colour", "black"))
    tcol = textcolor or textcolour or kwargs.get("color", kwargs.get("colour", "black"))
    shared_params = {
        "node": node,
        "label": label,
        "offset": offset,
        "offset_text": offset_text,
        "extend": extend,
        "align": align,
        "barsize": barsize,
        "barcolor": bcol,
        "textcolor": tcol,
        "angle": angle,
        "fontsize": fontsize,
        "horizontal": horizontal,
        "geom": geom,
        "hjust": hjust,
        **kwargs,
    }
    bar_layer = Layer("cladelab_bar", mapping=m, data=d, params=shared_params)
    txt_layer = Layer("cladelab_text", mapping=m, data=d, params=shared_params)
    return LayerList([bar_layer, txt_layer])


def geom_cladelabel(node: Any = None, label: Any = None, **kwargs: Any) -> LayerList:
    """Annotate a clade with a bar and label (R/geom_cladelabel.R)."""
    return geom_cladelab(node=node, label=label, **kwargs)


def geom_cladelabel2(node: Any = None, label: Any = None, **kwargs: Any) -> LayerList:
    """Annotate a clade in unrooted/circular layouts (R/geom_cladelabel2.R)."""
    return geom_cladelab(node=node, label=label, **kwargs)


def geom_strip(
    taxa1: Any = None,
    taxa2: Any = None,
    label: Any = None,
    offset: float = 0.0,
    offset_text: float = 0.0,
    align: bool = True,
    barsize: float = 0.5,
    extend: Any = 0.0,
    fontsize: float = 3.88,
    angle: Any = 0.0,
    geom: str = "text",
    hjust: float = 0.0,
    fill: Optional[str] = None,
    family: str = "sans",
    parse: bool = False,
    horizontal: bool = True,
    data: Any = None,
    mapping: Any = None,
    **kwargs: Any,
) -> LayerList:
    """Add strip bar + optional label spanning between taxa1 and taxa2 (R/geom_strip.R)."""
    offset_text = kwargs.pop("offset.text", offset_text)
    m, d = _extract_mapping_and_data(mapping, data)
    shared = {
        "taxa1": taxa1,
        "taxa2": taxa2,
        "label": label,
        "offset": offset,
        "offset_text": offset_text,
        "align": align,
        "barsize": barsize,
        "extend": extend,
        "fontsize": fontsize,
        "angle": angle,
        "geom": geom,
        "hjust": hjust,
        "fill": fill,
        "family": family,
        "parse": parse,
        "horizontal": horizontal,
        **kwargs,
    }
    bar_layer = Layer("strip_bar", mapping=m, data=d, params=shared)
    if label is not None or "label" in m:
        txt_layer = Layer("strip_text", mapping=m, data=d, params=shared)
        return LayerList([bar_layer, txt_layer])
    return LayerList([bar_layer])


def geom_striplab(
    taxa1: Any = None,
    taxa2: Any = None,
    label: Any = None,
    data: Any = None,
    mapping: Any = None,
    **kwargs: Any,
) -> LayerList:
    """Add strip bar and label via aesthetic mapping or direct taxa1/taxa2 (R/geom_strip.R)."""
    return geom_strip(taxa1=taxa1, taxa2=taxa2, label=label, data=data, mapping=mapping, **kwargs)


def geom_taxalink(
    taxa1: Any = None,
    taxa2: Any = None,
    offset: Optional[float] = None,
    outward: Union[str, bool] = "auto",
    curvature: float = 0.5,
    color: str = "black",
    linetype: Any = 1,
    size: float = 0.5,
    alpha: float = 1.0,
    arrow: Any = None,
    hratio: float = 1.0,
    ncp: int = 10,
    data: Any = None,
    mapping: Any = None,
    **kwargs: Any,
) -> Layer:
    """Link related taxa with curves (R/geom_taxalink.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    params = {
        "taxa1": taxa1,
        "taxa2": taxa2,
        "offset": offset,
        "outward": outward,
        "curvature": curvature,
        "color": color,
        "linetype": linetype,
        "size": size,
        "alpha": alpha,
        "arrow": arrow,
        "hratio": hratio,
        "ncp": ncp,
        **kwargs,
    }
    return Layer("taxalink", mapping=m, data=d, params=params)


def geom_tanglelink(
    assoc: Any = None,
    mapping: Any = None,
    data: Any = None,
    curvature: float = 0.0,
    color: str = "grey50",
    alpha: float = 0.8,
    linetype: Any = 1,
    linewidth: float = 0.5,
    size: Optional[float] = None,
    arrow: Any = None,
    **kwargs: Any,
) -> Layer:
    """Draw connecting links between two trees in a tanglegram / ggdoubletree (R/tanglegram.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    lw = float(size) if size is not None else float(linewidth)
    params = {
        "assoc": assoc,
        "curvature": curvature,
        "color": color,
        "alpha": alpha,
        "linetype": linetype,
        "linewidth": lw,
        "arrow": arrow,
        **kwargs,
    }
    return Layer("tanglelink", mapping=m, data=d, params=params)


def geom_treescale(
    x: Optional[float] = None,
    y: Optional[float] = None,
    width: Optional[float] = None,
    offset: Optional[float] = None,
    offset_label: Optional[float] = None,
    color: str = "black",
    linesize: float = 0.5,
    fontsize: float = 3.88,
    family: str = "sans",
    **kwargs: Any,
) -> Layer:
    """Add evolutionary distance scale bar (R/geom_treescale.R)."""
    offset_label = kwargs.pop("offset.label", offset_label)
    params = {
        "x": x,
        "y": y,
        "width": width,
        "offset": offset,
        "offset_label": offset_label,
        "color": color,
        "linesize": linesize,
        "fontsize": fontsize,
        "family": family,
        **kwargs,
    }
    return Layer("treescale", params=params)


def geom_balance(
    node: Union[int, str],
    fill: str = "steelblue",
    color: str = "white",
    size: float = 0.5,
    alpha: float = 0.5,
    extend: float = 0.0,
    extendto: Optional[float] = None,
    **kwargs: Any,
) -> Layer:
    """Highlight the two child subclades of an internal node (R/geom_balance.R)."""
    params = {
        "node": node,
        "fill": fill,
        "color": color,
        "size": size,
        "alpha": alpha,
        "extend": extend,
        "extendto": extendto,
        **kwargs,
    }
    return Layer("balance", params=params)


def geom_range(
    range: str = "length_0.95_HPD",
    center: str = "auto",
    color: str = "red",
    size: float = 2.0,
    alpha: float = 0.6,
    mapping: Any = None,
    data: Any = None,
    **kwargs: Any,
) -> Layer:
    """Add horizontal uncertainty/HPD bars on internal nodes (R/geom_range.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    params = {"range": range, "center": center, "color": color, "size": size, "alpha": alpha, **kwargs}
    return Layer("range", mapping=m, data=d, params=params)


def geom_rootedge(rootedge: Optional[float] = None, **kwargs: Any) -> Layer:
    """Add root edge line segment to a tree plot (R/geom_rootedge.R)."""
    return Layer("rootedge", params={"rootedge": rootedge, **kwargs})


def geom_zoom_clade(node: Union[int, str], xexpand: Optional[float] = None) -> PlotModifier:
    """Zoom into a specific clade (R/geom-zoom-clade.R)."""
    from .manipulate import zoomClade

    class _ZoomCladeMod(PlotModifier):
        def apply(self, plot: Any) -> Any:
            return zoomClade(plot, node=node, xexpand=xexpand)

    return _ZoomCladeMod()


def geom_motif(mapping: Any = None, data: Any = None, on: Optional[str] = None, label: Optional[str] = None, align: str = "left", **kwargs: Any) -> Layer:
    """Draw gene/domain motifs aligned with tree tips (R/experimental_function.R)."""
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("motif", mapping=m, data=d, params={"on": on, "label": label, "align": align, **kwargs})


def geom_inset(
    insets: Dict[Any, Any],
    width: float = 0.1,
    height: float = 0.1,
    hjust: float = 0.0,
    vjust: float = 0.0,
    x: str = "node",
    reverse_x: bool = False,
    reverse_y: bool = False,
    **kwargs: Any,
) -> Layer:
    """Add subplots (nodebar / nodepie) as insets on tree nodes (R/inset.R)."""
    params = {
        "insets": dict(insets),
        "width": width,
        "height": height,
        "hjust": hjust,
        "vjust": vjust,
        "x": x,
        "reverse_x": reverse_x,
        "reverse_y": reverse_y,
        **kwargs,
    }
    return Layer("inset", params=params)


# Standard ggplot2 geoms used in facet_plot / geom_facet / gheatmap / msaplot
def geom_bar(mapping: Any = None, data: Any = None, stat: str = "identity", **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("bar", mapping=m, data=d, params=kwargs, stat=stat)


def geom_col(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("col", mapping=m, data=d, params=kwargs)


def geom_tile(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("tile", mapping=m, data=d, params=kwargs)


def geom_raster(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    return geom_tile(mapping=mapping, data=data, **kwargs)


def geom_line(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("line", mapping=m, data=d, params=kwargs)


def geom_path(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("path", mapping=m, data=d, params=kwargs)


def geom_area(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("area", mapping=m, data=d, params=kwargs)


def geom_ribbon(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("ribbon", mapping=m, data=d, params=kwargs)


def geom_rect(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("rect", mapping=m, data=d, params=kwargs)


def geom_polygon(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("polygon", mapping=m, data=d, params=kwargs)


def geom_curve(mapping: Any = None, data: Any = None, curvature: float = 0.5, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("curve", mapping=m, data=d, params={"curvature": curvature, **kwargs})


def geom_vline(xintercept: Optional[float] = None, mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("vline", mapping=m, data=d, params={"xintercept": xintercept, **kwargs})


def geom_hline(yintercept: Optional[float] = None, mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("hline", mapping=m, data=d, params={"yintercept": yintercept, **kwargs})


def geom_abline(intercept: float = 0.0, slope: float = 1.0, **kwargs: Any) -> Layer:
    return Layer("abline", params={"intercept": intercept, "slope": slope, **kwargs})


def geom_boxplot(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("boxplot", mapping=m, data=d, params=kwargs)


def geom_violin(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("violin", mapping=m, data=d, params=kwargs)


def geom_density(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("density", mapping=m, data=d, params=kwargs)


def geom_histogram(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("histogram", mapping=m, data=d, params=kwargs)


def geom_errorbar(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("errorbar", mapping=m, data=d, params=kwargs)


def geom_errorbarh(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("errorbarh", mapping=m, data=d, params=kwargs)


def geom_linerange(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("linerange", mapping=m, data=d, params=kwargs)


def geom_pointrange(mapping: Any = None, data: Any = None, **kwargs: Any) -> Layer:
    m, d = _extract_mapping_and_data(mapping, data)
    return Layer("pointrange", mapping=m, data=d, params=kwargs)


# Interactive geom aliases and ggproto objects (from R/geom_interactive.R, R/AllClasses.R)
_INTERACTIVE_OPTIONS: Dict[str, Any] = {"enabled": False}


def ggtree_set_interactive(**kwargs: Any) -> Dict[str, Any]:
    """Enable interactive ggiraph tooltips/hover mode globally (R/geom_interactive.R)."""
    _INTERACTIVE_OPTIONS["enabled"] = True
    _INTERACTIVE_OPTIONS.update(kwargs)
    return dict(_INTERACTIVE_OPTIONS)


def ggtree_unset_interactive() -> Dict[str, Any]:
    """Disable interactive ggiraph mode globally (R/geom_interactive.R)."""
    _INTERACTIVE_OPTIONS["enabled"] = False
    return dict(_INTERACTIVE_OPTIONS)


class _GGProto:
    """Lightweight ggproto-compatible metadata & compute wrapper matching ggplot2/ggtree classes."""

    def __init__(self, name: str, required_aes: Sequence[str] = (), default_aes: Optional[Dict[str, Any]] = None, compute_fn: Any = None):
        self.__name__ = name
        self.required_aes = tuple(required_aes)
        self.default_aes = dict(default_aes or {})
        self._compute_fn = compute_fn

    def compute_panel(self, data: DataFrame, params: Optional[Dict[str, Any]] = None, **kwargs: Any) -> DataFrame:
        if self._compute_fn is not None:
            return self._compute_fn(data, params or {}, **kwargs)
        return data

    def __repr__(self) -> str:
        return f"<ggproto object: {self.__name__}>"


StatBalance = _GGProto("StatBalance", required_aes=("node",), compute_fn=lambda df, p, **kw: _compute_balance(df, p, **kw))
GeomHilightEncircle = _GGProto("GeomHilightEncircle", required_aes=("node",), default_aes={"fill": "steelblue", "alpha": 0.5})
GeomInteractiveCurvelink = _GGProto("GeomInteractiveCurvelink", required_aes=("x", "y", "xend", "yend"))
GeomInteractiveHilightEncircle = _GGProto("GeomInteractiveHilightEncircle", required_aes=("node",))
GeomInteractiveHilightRect = _GGProto("GeomInteractiveHilightRect", required_aes=("xmin", "xmax", "ymin", "ymax"))
GeomInteractiveLabelGGtree = _GGProto("GeomInteractiveLabelGGtree", required_aes=("x", "y", "label"))
GeomInteractivePointGGtree = _GGProto("GeomInteractivePointGGtree", required_aes=("x", "y"))
GeomInteractiveSegmentGGtree = _GGProto("GeomInteractiveSegmentGGtree", required_aes=("x", "y", "xend", "yend"))
GeomInteractiveTextGGtree = _GGProto("GeomInteractiveTextGGtree", required_aes=("x", "y", "label"))

geom_point_interactive = geom_point2
geom_point2_interactive = geom_point2
geom_text_interactive = geom_text2
geom_text2_interactive = geom_text2
geom_label_interactive = geom_label2
geom_label2_interactive = geom_label2
geom_segment_interactive = geom_segment2
geom_segment2_interactive = geom_segment2
geom_tiplab_interactive = geom_tiplab
geom_tiplab2_interactive = geom_tiplab2
geom_nodelab_interactive = geom_nodelab
geom_nodelab2_interactive = geom_nodelab2
geom_cladelab_interactive = geom_cladelab
geom_cladelabel_interactive = geom_cladelabel
geom_hilight_interactive = geom_hilight
geom_highlight_interactive = geom_hilight
geom_strip_interactive = geom_strip
geom_striplab_interactive = geom_striplab
geom_taxalink_interactive = geom_taxalink
geom_tree_interactive = geom_tree

