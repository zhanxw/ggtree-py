"""
Composite, multi-panel, and multi-tree functions for ggtree:
gheatmap, get_heatmap_column_position, msaplot, facet_plot, geom_facet,
facet_data, facet_labeller, facet_widths, inset, nodebar, nodepie,
ggdensitree, as_tanglegram, ggdoubletree, multiplot, plot_list.
"""

from __future__ import annotations

import copy
import math
import warnings
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union
import numpy as np

from .aes import Aes, aes
from .dataframe import (
    ColumnSeries,
    DataFrame,
    Factor,
    TblTree,
    _is_scalar_na,
    bind_rows,
    factor,
    left_join,
)
from .geoms import (
    Layer,
    geom_inset,
    geom_segment,
    geom_tanglelink,
    geom_text,
    geom_tile,
    geom_tiplab,
    geom_tree,
)
from .layout import (
    calculate_angle,
    calculate_branch_mid,
    fortify,
    getYcoord,
    re_assign_ycoord_df,
    tip_y_alignment,
)
from .themes_scales import (
    PlotModifier,
    scale_color_gradient,
    scale_fill_gradient,
    scale_fill_manual,
)
from .tidytree import child, getRoot, offspring
from .tree import MultiPhylo, Nnode2, Ntip, Phylo, TreeData, as_phylo, read_fasta


def _warn_misaligned_tips(p: Any, caller: str) -> None:
    """Port of warn_misaligned_tips in R/layout.R."""
    status = tip_y_alignment(p)
    if status == "duplicated":
        warnings.warn(
            f"{caller}(): some tips share the same vertical position (e.g. layout = 'tidy'); "
            "tip-aligned rows will overlap."
        )
    elif status in ("irregular", "uneven"):
        warnings.warn(
            f"{caller}(): tips are not evenly spaced on the y axis; "
            "fixed-height tiles may not align with tips."
        )


# 1. gheatmap & get_heatmap_column_position (R/gheatmap.R)
def gheatmap(
    p: Any,
    data: Any,
    offset: float = 0.0,
    width: float = 1.0,
    low: str = "green",
    high: str = "red",
    color: Optional[str] = "white",
    colnames: bool = True,
    colnames_position: str = "bottom",
    colnames_angle: float = 0.0,
    colnames_level: Optional[Sequence[str]] = None,
    colnames_offset_x: float = 0.0,
    colnames_offset_y: float = 0.0,
    font_size: float = 4.0,
    family: str = "",
    hjust: float = 0.5,
    legend_title: str = "value",
    custom_column_labels: Optional[Union[Sequence[str], Callable[[Sequence[str]], Sequence[str]]]] = None,
    **kwargs: Any,
) -> Any:
    """
    Append a heatmap of a matrix/DataFrame to the right side of a phylogenetic tree plot
    (exact port of ggtree::gheatmap in R/gheatmap.R).
    """
    font_size = kwargs.pop("font.size", font_size)
    if colnames_position not in ("bottom", "top"):
        raise ValueError("colnames_position should be 'bottom' or 'top'.")

    _warn_misaligned_tips(p, "gheatmap")
    out = p.copy()

    if isinstance(data, np.ndarray):
        arr = data
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        cols_auto = [f"V{j + 1}" for j in range(arr.shape[1])]
        df_mat = DataFrame({c: arr[:, j].tolist() for j, c in enumerate(cols_auto)})
    elif not isinstance(data, DataFrame):
        df_mat = DataFrame(data)
    else:
        df_mat = data.copy()

    if len(df_mat) == 0 or len(df_mat.columns) == 0:
        return out

    tree_df = out.data
    tips_df = tree_df[tree_df["isTip"]] if "isTip" in tree_df else tree_df
    tip_records = [
        (str(r["label"]), float(r["y"]))
        for r in tips_df.to_records()
        if not _is_scalar_na(r.get("label")) and not _is_scalar_na(r.get("y"))
    ]
    tip_y_map = dict(tip_records)

    # Match rows of df_mat to tree tips via rownames or first column if rownames don't match
    rnames = df_mat.rownames
    cols = list(df_mat.columns)
    if not any(rn in tip_y_map for rn in rnames):
        if len(df_mat) == len(tip_records) and isinstance(data, np.ndarray):
            # Unnamed numpy matrix matching number of tips -> align by tip label order
            rnames = [lbl for lbl, _ in tip_records]
        elif len(cols) > 1:
            first_col_vals = [str(v) for v in df_mat[cols[0]].tolist()]
            if any(v in tip_y_map for v in first_col_vals):
                rnames = first_col_vals
                cols = cols[1:]

    if colnames_level is not None:
        cols = [c for c in colnames_level if c in cols]

    rname_to_row = {rn: idx for idx, rn in enumerate(rnames)}
    missing_tips = [lbl for lbl, _ in tip_records if lbl not in rname_to_row]
    if missing_tips:
        warnings.warn(
            "The following tree labels are missing from `data` row names and will be filled with NA: "
            + ", ".join(missing_tips[:10])
        )

    xs = [float(v) for v in tree_df["x"].tolist() if not _is_scalar_na(v)]
    ys = [float(v) for v in tree_df["y"].tolist() if not _is_scalar_na(v)]
    min_x, max_x = (min(xs), max(xs)) if xs else (0.0, 1.0)
    span_x = max(max_x - min_x, 1e-6)

    ncol = max(len(cols), 1)
    col_width = float(width) * span_x / float(ncol)
    start_x = max_x + float(offset) + col_width / 2.0

    col_positions: List[float] = [start_x + j * col_width for j in range(len(cols))]

    tile_rows: List[Dict[str, Any]] = []
    recs = df_mat.to_records()
    for lbl, ty in tip_records:
        rec = recs[rname_to_row[lbl]] if lbl in rname_to_row else None
        for j, cname in enumerate(cols):
            val = rec.get(cname, None) if rec is not None else None
            tile_rows.append({
                "x": col_positions[j],
                "y": ty,
                "variable": cname,
                "value": val,
                "width": col_width,
                "height": 1.0,
                "row_label": lbl,
            })

    tile_df = DataFrame(tile_rows)
    tile_params = {"color": color, "width": col_width, "height": 1.0, "inherit.aes": False}
    out = out + geom_tile(
        mapping=aes(x="x", y="y", fill="value"),
        data=tile_df,
        **tile_params,
    )

    # Determine if values are numeric or categorical to add appropriate default fill scale
    vals_non_na = [r["value"] for r in tile_rows if not _is_scalar_na(r["value"])]
    is_num = len(vals_non_na) > 0 and all(isinstance(v, (int, float, np.number)) and not isinstance(v, bool) for v in vals_non_na)
    if is_num and low is not None and high is not None:
        out = out + scale_fill_gradient(low=low, high=high, name=legend_title)

    # Column labels
    min_y = min(ys) if ys else 1.0
    max_y = max(ys) if ys else 1.0
    y_pos = (min_y - 0.6) if colnames_position == "bottom" else (max_y + 0.6)
    y_pos += float(colnames_offset_y)

    display_labels = list(cols)
    if custom_column_labels is not None:
        if callable(custom_column_labels):
            display_labels = [str(x) for x in custom_column_labels(cols)]
        elif isinstance(custom_column_labels, dict):
            # Support both {orig_col: new_label} and R named vector c(new_label = "orig_col")
            if any(str(k) in cols for k in custom_column_labels.keys()):
                cmap = {str(k): str(v) for k, v in custom_column_labels.items()}
            else:
                cmap = {str(v): str(k) for k, v in custom_column_labels.items()}
            display_labels = [cmap.get(c, c) for c in cols]
        else:
            custom_list = [str(x) for x in custom_column_labels]
            if len(custom_list) <= len(cols):
                display_labels[: len(custom_list)] = custom_list

    col_pos_df = DataFrame({
        "x": [cx + float(colnames_offset_x) for cx in col_positions],
        "y": [y_pos] * len(cols),
        "label": display_labels,
    })
    data_axis_info = {
        "from": float(start_x - col_width / 2.0),
        "to": float(start_x + (len(cols) - 0.5) * col_width),
        "y": float(y_pos),
        "custom_labels": list(display_labels),
        "column_position": col_pos_df.copy(),
    }
    if hasattr(out, "attrs") and isinstance(out.attrs, dict):
        out.attrs["data_axis"] = data_axis_info
    out.data.attrs["column_position"] = col_pos_df

    if colnames:
        out = out + geom_text(
            mapping=aes(x="x", y="y", label="label"),
            data=col_pos_df,
            angle=colnames_angle,
            size=font_size,
            family=family,
            hjust=hjust,
            **{"inherit.aes": False},
        )

    return out


def get_heatmap_column_position(treeview: Any, by: str = "bottom") -> DataFrame:
    """Retrieve column label coordinates produced by gheatmap (R/experimental_function.R)."""
    pos = None
    if hasattr(treeview, "attrs") and isinstance(treeview.attrs, dict):
        da = treeview.attrs.get("data_axis")
        if isinstance(da, dict):
            pos = da.get("column_position")
    if pos is None:
        pos = treeview.data.attrs.get("column_position")
    if pos is None:
        return DataFrame(columns=["x", "y", "label"])
    out = pos.copy()
    ys = [float(v) for v in treeview.data["y"].tolist() if not _is_scalar_na(v)]
    if ys:
        y_val = (min(ys) - 0.6) if by == "bottom" else (max(ys) + 0.6)
        out["y"] = [y_val] * len(out)
    return out


# 2. msaplot (R/msaplot.R)
def msaplot(
    p: Any,
    fasta: Any,
    offset: float = 0.0,
    width: float = 1.0,
    color: Optional[Sequence[str]] = None,
    window: Optional[Tuple[int, int]] = None,
    bg_line: bool = True,
    height: float = 0.8,
) -> Any:
    """
    Visualize multiple sequence alignment alongside a phylogenetic tree
    (exact port of ggtree::msaplot in R/msaplot.R).
    """
    _warn_misaligned_tips(p, "msaplot")
    out = p.copy()

    if isinstance(fasta, str):
        seqs = read_fasta(fasta)
    elif isinstance(fasta, dict):
        seqs = {str(k): str(v) for k, v in fasta.items()}
    else:
        seqs = {}

    if not seqs:
        return out

    first_seq = next(iter(seqs.values()))
    seq_len_total = len(first_seq)
    if window is not None and len(window) >= 2:
        w_start = max(1, int(window[0]))
        w_end = min(seq_len_total, int(window[1]))
    else:
        w_start, w_end = 1, seq_len_total

    slice_len = max(w_end - w_start + 1, 1)
    tree_df = out.data
    tips_df = tree_df[tree_df["isTip"]] if "isTip" in tree_df else tree_df
    tip_y_map = {str(r["label"]): float(r["y"]) for r in tips_df.to_records() if not _is_scalar_na(r.get("label"))}

    xs = [float(v) for v in tree_df["x"].tolist() if not _is_scalar_na(v)]
    min_x, max_x = (min(xs), max(xs)) if xs else (0.0, 1.0)
    span_x = max(max_x - min_x, 1e-6)
    step_x = float(width) * span_x / float(slice_len)
    start_x = max_x + float(offset) + step_x / 2.0

    msa_rows: List[Dict[str, Any]] = []
    line_rows: List[Dict[str, Any]] = []
    for name, s_str in seqs.items():
        if name not in tip_y_map:
            continue
        ty = tip_y_map[name]
        sub_s = s_str[w_start - 1 : w_end]
        if bg_line:
            line_rows.append({
                "x": start_x,
                "xend": start_x + (len(sub_s) - 1) * step_x,
                "y": ty,
                "yend": ty,
                "name": name,
            })
        for pos_idx, ch in enumerate(sub_s):
            if ch in ("-", ".", "?", "X", "x"):
                continue
            msa_rows.append({
                "x": start_x + pos_idx * step_x,
                "y": ty,
                "site": w_start + pos_idx,
                "base": ch.upper(),
                "name": name,
                "width": step_x,
                "height": float(height),
            })

    if bg_line and line_rows:
        out = out + geom_segment(
            mapping=aes(x="x", xend="xend", y="y", yend="yend"),
            data=DataFrame(line_rows),
            color="grey70",
            size=0.3,
            inherit_aes=False,
        )
    if msa_rows:
        out = out + geom_tile(
            mapping=aes(x="x", y="y", fill="base"),
            data=DataFrame(msa_rows),
            width=step_x,
            height=float(height),
            inherit_aes=False,
        )
        if color is not None:
            out = out + scale_fill_manual(values=list(color))
    return out


# 3. facet_plot & geom_facet (R/facet_plot.R & R/facet-utilities.R)
class FacetLayer(PlotModifier):
    """Modifier returned by geom_facet(...) that aligns external data with tree tips in a named panel."""

    def __init__(
        self,
        panel: str,
        data: Any,
        geom: Callable[..., Any],
        mapping: Optional[Aes] = None,
        **kwargs: Any,
    ):
        self.panel = str(panel)
        self.data = data
        self.geom = geom
        self.mapping = mapping.copy() if isinstance(mapping, Aes) else Aes(**(mapping or {}))
        self.params = dict(kwargs)

    def apply(self, plot: Any) -> Any:
        _warn_misaligned_tips(plot, "facet_plot")
        out = plot.copy()
        if self.panel not in out.panels:
            out.panels.append(self.panel)

        df_ext = self.data.copy() if isinstance(self.data, DataFrame) else DataFrame(self.data)
        if len(df_ext) > 0:
            # Join with tip y coordinates from plot.data
            tree_df = out.data
            tips_df = tree_df[tree_df["isTip"]] if "isTip" in tree_df else tree_df
            tip_y_map = {str(r["label"]): float(r["y"]) for r in tips_df.to_records() if not _is_scalar_na(r.get("label"))}
            id_col = "label" if "label" in df_ext.columns else df_ext.columns[0]
            ys_aligned = [tip_y_map.get(str(v), None) for v in df_ext[id_col].tolist()]
            df_ext = df_ext.copy()
            if "y" not in df_ext:
                df_ext["y"] = ys_aligned
            df_ext[".panel"] = [self.panel] * len(df_ext)

        m = self.mapping.copy()
        if "y" not in m:
            m["y"] = "y"

        layer_obj = self.geom(mapping=m, data=df_ext, panel=self.panel, **self.params)
        if isinstance(layer_obj, Layer):
            layer_obj.params["panel"] = self.panel
        out.facet_data_store[self.panel] = df_ext
        out = out + layer_obj
        return out


def geom_facet(
    panel: Optional[str] = None,
    data: Any = None,
    geom: Callable[..., Any] = geom_tile,
    mapping: Optional[Aes] = None,
    **kwargs: Any,
) -> FacetLayer:
    """Add an associated data panel aligned to the tree tips (R/facet_plot.R)."""
    if panel is None:
        panel = "Trait"
    return FacetLayer(panel=panel, data=data, geom=geom, mapping=mapping, **kwargs)


def facet_plot(
    p: Any,
    panel: Optional[str] = None,
    data: Any = None,
    geom: Callable[..., Any] = geom_tile,
    mapping: Optional[Aes] = None,
    **kwargs: Any,
) -> Any:
    """Plot associated data in an adjacent facet panel aligned to tree tips (R/facet_plot.R)."""
    return p + geom_facet(panel=panel, data=data, geom=geom, mapping=mapping, **kwargs)


def facet_data(p: Any, panel: str) -> DataFrame:
    """Extract the data frame associated with a specific facet panel (R/facet-utilities.R)."""
    if panel == "Tree":
        return p.data
    return p.facet_data_store.get(panel, DataFrame())


def facet_labeller(p: Any, label: Dict[str, str]) -> Any:
    """Rename facet panel header labels (R/facet-utilities.R)."""
    out = p.copy()
    out.panel_labels.update(label)
    return out


def facet_widths(p: Any, widths: Union[Sequence[float], Dict[str, float]]) -> Any:
    """Set relative widths of facet panels (R/facet-utilities.R)."""
    out = p.copy()
    if isinstance(widths, dict):
        out.panel_widths.update(widths)
    else:
        for pan, w in zip(out.panels, widths):
            out.panel_widths[pan] = float(w)
    return out


# 4. inset, nodebar, nodepie (R/inset.R)
def nodebar(
    data: DataFrame,
    cols: Union[Sequence[int], Sequence[str]],
    color: Optional[Dict[str, str]] = None,
    alpha: float = 1.0,
    position: str = "stack",
) -> Dict[int, Dict[str, Any]]:
    """
    Generate bar chart inset definitions for tree nodes (exact port of ggtree::nodebar in R/inset.R).
    """
    df = data if isinstance(data, DataFrame) else DataFrame(data)
    if isinstance(cols, slice):
        col_names = list(df.columns[cols])
    else:
        col_names = [df.columns[int(c) - 1] if isinstance(c, int) else str(c) for c in cols]
    insets: Dict[int, Dict[str, Any]] = {}
    for r in df.to_records():
        nid = int(r.get("node", 0))
        vals = {c: float(r.get(c, 0.0) or 0.0) for c in col_names}
        insets[nid] = {
            "type": "bar",
            "node": nid,
            "values": vals,
            "color": color,
            "alpha": alpha,
            "position": position,
        }
    return insets


def nodepie(
    data: DataFrame,
    cols: Union[Sequence[int], Sequence[str]],
    color: Optional[Dict[str, str]] = None,
    alpha: float = 1.0,
) -> Dict[int, Dict[str, Any]]:
    """
    Generate pie chart inset definitions for tree nodes (exact port of ggtree::nodepie in R/inset.R).
    """
    df = data if isinstance(data, DataFrame) else DataFrame(data)
    if isinstance(cols, slice):
        col_names = list(df.columns[cols])
    else:
        col_names = [df.columns[int(c) - 1] if isinstance(c, int) else str(c) for c in cols]
    insets: Dict[int, Dict[str, Any]] = {}
    for r in df.to_records():
        nid = int(r.get("node", 0))
        vals = {c: float(r.get(c, 0.0) or 0.0) for c in col_names}
        insets[nid] = {
            "type": "pie",
            "node": nid,
            "values": vals,
            "color": color,
            "alpha": alpha,
        }
    return insets


def inset(
    tree_view: Any,
    insets: Dict[Any, Any],
    width: float = 0.1,
    height: float = 0.1,
    hjust: float = 0.0,
    vjust: float = 0.0,
    x: str = "node",
    reverse_x: bool = False,
    reverse_y: bool = False,
) -> Any:
    """Add nodebar / nodepie insets to a tree plot (R/inset.R)."""
    return tree_view + geom_inset(
        insets=insets,
        width=width,
        height=height,
        hjust=hjust,
        vjust=vjust,
        x=x,
        reverse_x=reverse_x,
        reverse_y=reverse_y,
    )


# 5. ggdensitree (R/ggdensitree.R)
def ggdensitree(
    data: Any = None,
    mapping: Optional[Aes] = None,
    layout: str = "slanted",
    tip_order: Union[str, Sequence[str]] = "mode",
    align_tips: bool = True,
    jitter: float = 0.0,
    **kwargs: Any,
) -> Any:
    """
    Draw multiple overlapping trees (densitree) with aligned tip labels
    (exact port of ggtree::ggdensitree in R/ggdensitree.R).
    """
    from .plot import GGTree

    tip_order = kwargs.pop("tip.order", tip_order)
    align_tips = kwargs.pop("align.tips", align_tips)

    if data is None:
        return GGTree(layout=layout, mapping=mapping)

    if isinstance(data, DataFrame):
        trees_f = [data.copy()]
    elif isinstance(data, (list, tuple, MultiPhylo)):
        trees_f = []
        for idx, item in enumerate(data, 1):
            if isinstance(item, DataFrame):
                d = item.copy()
            else:
                d = fortify(item, layout=layout)
            if "tree" not in d:
                d["tree"] = str(idx)
            trees_f.append(d)
    else:
        trees_f = [fortify(data, layout=layout)]

    if not trees_f:
        return GGTree(layout=layout, mapping=mapping)

    # Determine consensus / explicit tip_order
    if isinstance(tip_order, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(tip_order, (str, bytes)):
        ordered_labels = [str(x) for x in tip_order]
    else:
        first_tips = trees_f[0][trees_f[0]["isTip"]].sort_values("y")
        ordered_labels = [str(v) for v in first_tips["label"].tolist() if not _is_scalar_na(v)]

    label_to_y = {lbl: float(i) for i, lbl in enumerate(ordered_labels, 1)}

    aligned_dfs: List[DataFrame] = []
    for d in trees_f:
        df_c = d.copy()
        if align_tips:
            xs = [float(v) for v in df_c["x"].tolist() if not _is_scalar_na(v)]
            max_x = max(xs) if xs else 0.0
            df_c["x"] = [float(v) - max_x for v in df_c["x"].tolist()]
            if "branch" in df_c:
                df_c["branch"] = [float(v) - max_x for v in df_c["branch"].tolist()]

        # Recompute y coordinates from shared tip_order
        nodes_col = [int(n) for n in df_c["node"].tolist()]
        istip_col = [bool(t) for t in df_c["isTip"].tolist()]
        labels_col = df_c["label"].tolist()
        new_y: List[Optional[float]] = []
        tip_nodes: List[int] = []
        for n, tflag, lbl in zip(nodes_col, istip_col, labels_col):
            if tflag:
                yv = label_to_y.get(str(lbl), 1.0)
                if jitter > 0:
                    yv += float(np.random.normal(0.0, float(jitter)))
                new_y.append(yv)
                tip_nodes.append(n)
            else:
                new_y.append(None)
        df_c["y"] = new_y
        df_c = re_assign_ycoord_df(df_c, tip_nodes)
        df_c = calculate_angle(df_c)
        aligned_dfs.append(df_c)

    combined = bind_rows(*aligned_dfs)
    combined.attrs["layout"] = layout
    p = GGTree(data=combined, mapping=mapping, layout=layout)
    tree_params = {k: v for k, v in kwargs.items() if k in ("alpha", "color", "colour", "size", "linewidth", "linetype")}
    p = p + geom_tree(layout=layout, **tree_params)
    return p


# 6. Tanglegram: as_tanglegram & ggdoubletree (R/tanglegram.R)
class _DiagDict(dict):
    """Dict with attribute access for tangle_optimize diagnostics."""

    def __getattr__(self, name: str) -> Any:
        if name in self:
            return self[name]
        raise AttributeError(name)


def _tip_y_map(df: DataFrame) -> Dict[str, float]:
    tips = df[df["isTip"]] if "isTip" in df else df
    return {
        str(r["label"]): float(r["y"])
        for r in tips.to_records()
        if not _is_scalar_na(r.get("label")) and not _is_scalar_na(r.get("y"))
    }


def _count_tangle_crossings(
    left_y: Dict[str, float],
    right_y: Dict[str, float],
    pairs: Sequence[Tuple[str, str]],
) -> int:
    pts = [(left_y[l], right_y[r]) for l, r in pairs if l in left_y and r in right_y]
    m = len(pts)
    crossings = 0
    for i in range(m):
        y1i, y2i = pts[i]
        for j in range(i + 1, m):
            dy1 = y1i - pts[j][0]
            dy2 = y2i - pts[j][1]
            if dy1 * dy2 < 0:
                crossings += 1
    return crossings


def _apply_child_order_to_df(df: DataFrame, root_id: int, kids: Dict[int, List[int]]) -> DataFrame:
    out = df.copy()
    nodes = [int(n) for n in out["node"].tolist()]
    is_tip = [bool(t) for t in out["isTip"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes)}

    ordered_tips: List[int] = []
    pre_stack = [root_id]
    while pre_stack:
        u = pre_stack.pop()
        uidx = id_to_idx[u]
        if is_tip[uidx]:
            ordered_tips.append(u)
        else:
            for ch in reversed(kids.get(u, [])):
                pre_stack.append(ch)

    new_y: List[Optional[float]] = [None] * len(out)
    for rank_i, tid in enumerate(ordered_tips, 1):
        new_y[id_to_idx[tid]] = float(rank_i)
    out["y"] = new_y
    out = re_assign_ycoord_df(out, ordered_tips)
    out = calculate_angle(out)
    return out


def _untangle_tree_against(
    df_target: DataFrame,
    fixed_tip_y: Dict[str, float],
    target_is_left: bool,
    pairs: Sequence[Tuple[str, str]],
    max_iter: int = 8,
) -> Tuple[DataFrame, int]:
    """
    Non-regressive untangling of df_target against fixed_tip_y:
    combines barycenter sweep + internal-node child permutation hill-climbing,
    only accepting moves that strictly reduce crossings.
    """
    best_df = df_target.copy()
    cur_target_y = _tip_y_map(best_df)
    best_crossings = (
        _count_tangle_crossings(cur_target_y, fixed_tip_y, pairs)
        if target_is_left
        else _count_tangle_crossings(fixed_tip_y, cur_target_y, pairs)
    )
    if best_crossings == 0:
        return best_df, 0

    to_fixed: Dict[str, List[str]] = {}
    for l_lbl, r_lbl in pairs:
        if target_is_left:
            to_fixed.setdefault(l_lbl, []).append(r_lbl)
        else:
            to_fixed.setdefault(r_lbl, []).append(l_lbl)

    target_scores = {
        t_lbl: sum(fixed_tip_y[f] for f in f_list if f in fixed_tip_y) / len([f for f in f_list if f in fixed_tip_y])
        for t_lbl, f_list in to_fixed.items()
        if any(f in fixed_tip_y for f in f_list)
    }

    iters = 0
    if target_scores:
        cand = _untangle_single_sweep(best_df, target_scores)
        cand_y = _tip_y_map(cand)
        cand_c = (
            _count_tangle_crossings(cand_y, fixed_tip_y, pairs)
            if target_is_left
            else _count_tangle_crossings(fixed_tip_y, cand_y, pairs)
        )
        iters += 1
        if cand_c < best_crossings:
            best_df = cand
            best_crossings = cand_c
            cur_target_y = cand_y

    if best_crossings == 0:
        return best_df, iters

    # Build internal tree structure ordered by current tip y
    nodes = [int(n) for n in best_df["node"].tolist()]
    parents = [int(p) for p in best_df["parent"].tolist()]
    is_tip = [bool(t) for t in best_df["isTip"].tolist()]
    labels = best_df["label"].tolist()
    ys_cur = [float(v) if not _is_scalar_na(v) else 0.0 for v in best_df["y"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes)}
    root_id = getRoot(best_df)

    kids: Dict[int, List[int]] = {}
    for p, n in zip(parents, nodes):
        if p != n:
            kids.setdefault(p, []).append(n)
    for u in kids:
        kids[u].sort(key=lambda c: ys_cur[id_to_idx[c]])

    # Local search over internal node child orders (guaranteed non-regressive)
    import itertools

    def _current_tip_y_map(kmap: Dict[int, List[int]]) -> Dict[str, float]:
        res: Dict[str, float] = {}
        rank = 1
        stk = [root_id]
        while stk:
            u = stk.pop()
            uidx = id_to_idx[u]
            if is_tip[uidx]:
                lbl = labels[uidx]
                if not _is_scalar_na(lbl):
                    res[str(lbl)] = float(rank)
                rank += 1
            else:
                for ch in reversed(kmap.get(u, [])):
                    stk.append(ch)
        return res

    internal_nodes = [n for n in kids if len(kids[n]) >= 2]
    improved = True
    while improved and best_crossings > 0 and iters < max_iter:
        improved = False
        iters += 1
        for u in internal_nodes:
            ch_list = list(kids[u])
            if len(ch_list) <= 3:
                perms = list(itertools.permutations(ch_list))
            else:
                perms = [tuple(reversed(ch_list))] + [
                    tuple(ch_list[:k] + [ch_list[k + 1], ch_list[k]] + ch_list[k + 2:])
                    for k in range(len(ch_list) - 1)
                ]
            best_perm = tuple(ch_list)
            for perm in perms:
                if perm == tuple(ch_list):
                    continue
                kids[u] = list(perm)
                ty = _current_tip_y_map(kids)
                c_val = (
                    _count_tangle_crossings(ty, fixed_tip_y, pairs)
                    if target_is_left
                    else _count_tangle_crossings(fixed_tip_y, ty, pairs)
                )
                if c_val < best_crossings:
                    best_crossings = c_val
                    best_perm = perm
                    improved = True
            kids[u] = list(best_perm)
            if best_crossings == 0:
                break

    if improved or iters > 1:
        best_df = _apply_child_order_to_df(best_df, root_id, kids)
    return best_df, iters


def _untangle_single_sweep(df: DataFrame, target_tip_y: Dict[str, float]) -> DataFrame:
    """
    One post-order sweep that reorders children of internal nodes by mean target y
    of their descendant tips, then recomputes y coordinates (exact port of R/tanglegram.R).
    """
    out = df.copy()
    nodes = [int(n) for n in out["node"].tolist()]
    parents = [int(p) for p in out["parent"].tolist()]
    is_tip = [bool(t) for t in out["isTip"].tolist()]
    labels = out["label"].tolist()

    kids: Dict[int, List[int]] = {}
    for p, n in zip(parents, nodes):
        if p != n:
            kids.setdefault(p, []).append(n)

    id_to_idx = {n: i for i, n in enumerate(nodes)}
    root_id = getRoot(out)

    post_order: List[int] = []
    stack = [root_id]
    while stack:
        u = stack.pop()
        post_order.append(u)
        for ch in kids.get(u, []):
            stack.append(ch)
    post_order.reverse()

    node_score: Dict[int, float] = {}
    for u in post_order:
        uidx = id_to_idx[u]
        if is_tip[uidx]:
            lbl = str(labels[uidx]) if not _is_scalar_na(labels[uidx]) else ""
            if lbl in target_tip_y:
                node_score[u] = target_tip_y[lbl]
            else:
                node_score[u] = float(out["y"][uidx] or 0.0)
        else:
            chs = kids.get(u, [])
            scores = [node_score[c] for c in chs if c in node_score]
            node_score[u] = sum(scores) / len(scores) if scores else 0.0
            kids[u] = sorted(chs, key=lambda c: node_score.get(c, 0.0))

    return _apply_child_order_to_df(out, root_id, kids)


def _extract_tangle_input(obj: Any) -> Tuple[Any, List[Any]]:
    """Extract underlying tree/data object and any non-tree layers from GGTree/Phylo/TreeData."""
    if hasattr(obj, "layers") and hasattr(obj, "data"):
        non_tree = [
            lyr.copy()
            for lyr in obj.layers
            if getattr(lyr, "geom_type", "") not in ("tree", "tree_horizontal", "tree_vertical")
        ]
        return (obj.data.copy(), non_tree)
    return (obj, [])


class Tanglegram:
    """
    R 'gtangle' object produced by as_tanglegram(tr1, tr2, assoc=..., ...),
    fortifiable via fortify(tg) and plottable via ggtree(tg) + geom_tanglelink().
    """

    is_gtangle: bool = True

    def __init__(
        self,
        left: Any,
        right: Any,
        assoc: Optional[DataFrame] = None,
        gap: float = 0.08,
        optimize: bool = False,
        optimize_side: str = "right",
        mirror: bool = True,
        max_iter: int = 8,
    ):
        if optimize_side not in ("right", "left", "both"):
            raise ValueError("optimize_side must be one of 'right', 'left', or 'both'.")
        self.left, self.left_layers = _extract_tangle_input(left)
        self.right, self.right_layers = _extract_tangle_input(right)
        if assoc is not None and not isinstance(assoc, DataFrame):
            assoc = DataFrame(assoc)
        if assoc is not None and len(assoc.columns) < 2:
            raise ValueError("'assoc' must have at least 2 columns.")
        self.assoc = assoc
        self.gap = float(gap)
        self.optimize = bool(optimize)
        self.optimize_side = optimize_side
        self.mirror = bool(mirror)
        self.max_iter = int(max_iter)
        self.attrs: Dict[str, Any] = {
            "tangle_assoc": self.assoc,
            "tangle_gap": self.gap,
            "tangle_optimize": self.optimize,
            "tangle_optimize_side": self.optimize_side,
        }

    def copy(self) -> "Tanglegram":
        return Tanglegram(
            left=self.left,
            right=self.right,
            assoc=self.assoc.copy() if self.assoc is not None else None,
            gap=self.gap,
            optimize=self.optimize,
            optimize_side=self.optimize_side,
            mirror=self.mirror,
            max_iter=self.max_iter,
        )

    def fortify(
        self,
        layout: str = "rectangular",
        ladderize: bool = True,
        right: bool = False,
        branch_length: str = "branch.length",
        mrsd: Any = None,
        as_Date: bool = False,
        yscale: str = "none",
        **kwargs: Any,
    ) -> TblTree:
        opt_flag = bool(kwargs.pop("optimize", self.optimize))
        opt_side = str(kwargs.pop("optimize_side", self.optimize_side))
        if opt_side not in ("right", "left", "both"):
            raise ValueError("optimize_side must be one of 'right', 'left', or 'both'.")
        gap_val = float(kwargs.pop("gap", self.gap))
        mirror_val = bool(kwargs.pop("mirror", self.mirror))

        df1 = fortify(
            self.left,
            layout=layout,
            ladderize=ladderize,
            right=right,
            branch_length=branch_length,
            mrsd=mrsd,
            as_Date=as_Date,
            yscale=yscale,
            **kwargs,
        )
        df2 = fortify(
            self.right,
            layout=layout,
            ladderize=ladderize,
            right=right,
            branch_length=branch_length,
            mrsd=mrsd,
            as_Date=as_Date,
            yscale=yscale,
            **kwargs,
        )

        l_tip_y = _tip_y_map(df1)
        r_tip_y = _tip_y_map(df2)

        if self.assoc is not None and len(self.assoc) > 0:
            assoc_df = self.assoc.copy()
        else:
            l_ordered = [
                str(r["label"])
                for r in sorted(df1[df1["isTip"]].to_records(), key=lambda r: float(r["y"]))
                if not _is_scalar_na(r.get("label")) and str(r["label"]) in r_tip_y
            ]
            assoc_df = DataFrame({"left": l_ordered, "right": list(l_ordered)})

        c0, c1 = assoc_df.columns[0], assoc_df.columns[1]
        pairs: List[Tuple[str, str]] = []
        for rec in assoc_df.to_records():
            lv = "" if _is_scalar_na(rec.get(c0)) else str(rec.get(c0))
            rv = "" if _is_scalar_na(rec.get(c1)) else str(rec.get(c1))
            if lv and rv:
                pairs.append((lv, rv))

        before_c = _count_tangle_crossings(l_tip_y, r_tip_y, pairs)
        iters_done = 0

        if opt_flag and pairs:
            if opt_side in ("right", "both"):
                df2, n_it = _untangle_tree_against(df2, _tip_y_map(df1), target_is_left=False, pairs=pairs, max_iter=self.max_iter)
                iters_done += n_it
            if opt_side in ("left", "both"):
                df1, n_it = _untangle_tree_against(df1, _tip_y_map(df2), target_is_left=True, pairs=pairs, max_iter=self.max_iter)
                iters_done += n_it
            if opt_side == "both" and _count_tangle_crossings(_tip_y_map(df1), _tip_y_map(df2), pairs) > 0:
                df2, n_it = _untangle_tree_against(df2, _tip_y_map(df1), target_is_left=False, pairs=pairs, max_iter=self.max_iter)
                iters_done += n_it

        after_c = _count_tangle_crossings(_tip_y_map(df1), _tip_y_map(df2), pairs)
        diag = _DiagDict({
            "optimized": opt_flag,
            "optimize_side": opt_side,
            "iterations": iters_done,
            "before": before_c,
            "after": after_c,
        })

        # Offset and mirror right tree horizontally
        xs1 = [float(v) for v in df1["x"].tolist() if not _is_scalar_na(v)]
        xs2 = [float(v) for v in df2["x"].tolist() if not _is_scalar_na(v)]
        min_x1, max_x1 = (min(xs1), max(xs1)) if xs1 else (0.0, 1.0)
        min_x2, max_x2 = (min(xs2), max(xs2)) if xs2 else (0.0, 1.0)
        span1 = max(max_x1 - min_x1, 1e-6)
        span2 = max(max_x2 - min_x2, 1e-6)
        gap_dist = gap_val * (span1 + span2)

        df1 = df1.copy()
        df2 = df2.copy()
        df1["orig_node"] = list(df1["node"].tolist())
        df1["side"] = ["left"] * len(df1)
        df1["tree_id"] = ["left"] * len(df1)
        df1["tree"] = ["left"] * len(df1)

        node_offset = max(int(n) for n in df1["node"].tolist())
        df2["orig_node"] = list(df2["node"].tolist())
        df2["node"] = [int(n) + node_offset for n in df2["node"].tolist()]
        df2["parent"] = [int(p) + node_offset for p in df2["parent"].tolist()]

        def _transform_x2(val: float) -> float:
            if mirror_val:
                return max_x1 + gap_dist + (max_x2 - val)
            return max_x1 + gap_dist + (val - min_x2)

        df2["x"] = [_transform_x2(float(v)) if not _is_scalar_na(v) else None for v in df2["x"].tolist()]
        if "branch" in df2:
            df2["branch"] = [_transform_x2(float(v)) if not _is_scalar_na(v) else None for v in df2["branch"].tolist()]
        df2["side"] = ["right"] * len(df2)
        df2["tree_id"] = ["right"] * len(df2)
        df2["tree"] = ["right"] * len(df2)

        combined = bind_rows(df1, df2)

        # Build resolved tangle_assoc with node IDs and coordinates
        l_tips_map = {
            str(r["label"]): r
            for r in df1[df1["isTip"]].to_records()
            if not _is_scalar_na(r.get("label"))
        }
        r_tips_map = {
            str(r["label"]): r
            for r in df2[df2["isTip"]].to_records()
            if not _is_scalar_na(r.get("label"))
        }
        extra_cols = assoc_df.columns[2:]
        resolved_rows: List[Dict[str, Any]] = []
        for rec in assoc_df.to_records():
            lv = "" if _is_scalar_na(rec.get(c0)) else str(rec.get(c0))
            rv = "" if _is_scalar_na(rec.get(c1)) else str(rec.get(c1))
            if lv in l_tips_map and rv in r_tips_map:
                lr = l_tips_map[lv]
                rr = r_tips_map[rv]
                item: Dict[str, Any] = {
                    "left": lv,
                    "right": rv,
                    "left_node": int(lr["node"]),
                    "right_node": int(rr["node"]),
                    "x": float(lr["x"]),
                    "y": float(lr["y"]),
                    "xend": float(rr["x"]),
                    "yend": float(rr["y"]),
                }
                for ec in extra_cols:
                    item[ec] = rec.get(ec)
                resolved_rows.append(item)
        resolved_assoc = (
            DataFrame(resolved_rows)
            if resolved_rows
            else DataFrame(columns=["left", "right", "left_node", "right_node", "x", "y", "xend", "yend"])
        )

        combined.attrs["layout"] = "tanglegram"
        combined.attrs["tangle_assoc"] = resolved_assoc
        combined.attrs["tangle_optimize"] = diag
        combined.attrs["tangle_gap"] = gap_val
        combined.attrs["tangle_optimize_side"] = opt_side
        combined.attrs["tangle_transform"] = {
            "max_x1": max_x1,
            "max_x2": max_x2,
            "min_x2": min_x2,
            "gap_dist": gap_dist,
            "mirror": mirror_val,
            "node_offset": node_offset,
        }
        return TblTree(combined, attrs=combined.attrs)


def as_tanglegram(
    x: Any,
    y: Any = None,
    assoc: Any = None,
    gap: float = 0.08,
    optimize: bool = False,
    optimize_side: str = "right",
    mirror: bool = True,
    **kwargs: Any,
) -> Tanglegram:
    """
    Construct a Tanglegram ('gtangle') object from two trees (or GGTree plots / 2-tree list / cophylo)
    and an optional association table (exact port of ggtree::as_tanglegram in R/tanglegram.R).
    """
    if isinstance(x, Tanglegram) and y is None:
        out = x.copy()
        if assoc is not None:
            out.assoc = assoc if isinstance(assoc, DataFrame) else DataFrame(assoc)
        return out

    if y is None:
        if isinstance(x, dict) and "trees" in x and len(x["trees"]) >= 2:
            tr1, tr2 = x["trees"][0], x["trees"][1]
            if assoc is None and "assoc" in x:
                assoc = x["assoc"]
        elif hasattr(x, "trees") and len(getattr(x, "trees")) >= 2:
            tr1, tr2 = x.trees[0], x.trees[1]
            if assoc is None and hasattr(x, "assoc"):
                assoc = x.assoc
        elif isinstance(x, (list, tuple, MultiPhylo, TreeDataList)) and len(x) >= 2:
            tr1, tr2 = x[0], x[1]
        else:
            raise ValueError("as_tanglegram requires two trees (either (tr1, tr2) or a list of 2 trees).")
    else:
        tr1, tr2 = x, y
    return Tanglegram(
        left=tr1,
        right=tr2,
        assoc=assoc,
        gap=gap,
        optimize=optimize,
        optimize_side=optimize_side,
        mirror=mirror,
        **kwargs,
    )


as_tanglegram_cophylo = as_tanglegram
as_tanglegram_default = as_tanglegram
as_tanglegram_tanglegram = as_tanglegram
as_tanglegram.cophylo = as_tanglegram  # type: ignore[attr-defined]
as_tanglegram.default = as_tanglegram  # type: ignore[attr-defined]
as_tanglegram.tanglegram = as_tanglegram  # type: ignore[attr-defined]


def _replay_tangle_layer(lyr: Any, side_tag: str, tfm: Dict[str, Any]) -> Any:
    """Clone and adapt a non-tree layer from input GGTree to a tanglegram plot side."""
    from .geoms import Layer
    from .aes import Aes, SymExpr, eval_subset
    cloned = lyr.copy()
    m = cloned.mapping.copy() if isinstance(cloned.mapping, Aes) else Aes()
    saved_env = dict(getattr(m, "_env", {}))
    orig_sub = m.get("subset")
    side_expr = (
        SymExpr(
            f"(side == '{side_tag}' & ({orig_sub}))",
            lambda df, _s=orig_sub, _e=saved_env, _st=side_tag: (df["side"] == _st) & eval_subset(df, _s, extra_env=_e),
        )
        if orig_sub is not None
        else SymExpr(f"side == '{side_tag}'", lambda df, _st=side_tag: df["side"] == _st)
    )
    if cloned.data is None:
        m["subset"] = side_expr
        m._env = saved_env
        cloned.mapping = m
    if side_tag == "right":
        node_offset = int(tfm.get("node_offset", 0))
        if "node" in cloned.params and cloned.params["node"] is not None:
            nd = cloned.params["node"]
            if isinstance(nd, (int, np.integer)):
                cloned.params["node"] = int(nd) + node_offset
            elif isinstance(nd, (list, tuple, np.ndarray)):
                cloned.params["node"] = [int(v) + node_offset if isinstance(v, (int, np.integer)) else v for v in nd]
        if getattr(cloned, "geom_type", "") == "tiplab" and tfm.get("mirror", True):
            cloned.params["hjust"] = 1.0 - float(cloned.params.get("hjust", 0.0))
    return cloned


def ggdoubletree(
    tr1: Any,
    tr2: Any = None,
    assoc: Any = None,
    layout: str = "rectangular",
    gap: float = 0.08,
    optimize: bool = False,
    optimize_side: str = "right",
    mirror: bool = True,
    link_color: str = "grey50",
    link_alpha: float = 0.8,
    link_linetype: Any = 1,
    link_linewidth: Optional[float] = None,
    link_size: float = 0.5,
    tiplab: bool = True,
    **kwargs: Any,
) -> Any:
    """
    Create a side-by-side co-phylogeny / tanglegram plot of two trees or two GGTree plots
    (exact port of ggtree::ggdoubletree in R/tanglegram.R).
    """
    from .plot import ggtree
    tg = as_tanglegram(
        tr1,
        y=tr2,
        assoc=assoc,
        gap=gap,
        optimize=optimize,
        optimize_side=optimize_side,
        mirror=mirror,
    )
    df = fortify(tg, layout=layout)
    p = ggtree(df, **kwargs)
    tfm = df.attrs.get("tangle_transform", {})

    replayed_left_tiplab = False
    replayed_right_tiplab = False
    for lyr in tg.left_layers:
        p = p + _replay_tangle_layer(lyr, "left", tfm)
        if getattr(lyr, "geom_type", "") == "tiplab":
            replayed_left_tiplab = True
    for lyr in tg.right_layers:
        p = p + _replay_tangle_layer(lyr, "right", tfm)
        if getattr(lyr, "geom_type", "") == "tiplab":
            replayed_right_tiplab = True

    if tiplab:
        if not replayed_left_tiplab:
            p = p + geom_tiplab(aes(subset="side == 'left'"), hjust=0)
        if not replayed_right_tiplab:
            p = p + geom_tiplab(aes(subset="side == 'right'"), hjust=1 if mirror else 0)

    lw = float(link_linewidth) if link_linewidth is not None else float(link_size)
    if len(df.attrs.get("tangle_assoc", DataFrame())) > 0 or assoc is not None:
        p = p + geom_tanglelink(
            assoc=df.attrs.get("tangle_assoc"),
            color=link_color,
            alpha=link_alpha,
            linetype=link_linetype,
            linewidth=lw,
        )
    return p


# 7. multiplot & plot_list (R/multiplot.R & aplot::plot_list)
class PlotList:
    """Grid of multiple GGTree / ggplot objects (compatible with multiplot and aplot::plot_list)."""

    def __init__(
        self,
        plots: Sequence[Any],
        ncol: Optional[int] = None,
        nrow: Optional[int] = None,
        widths: Optional[Sequence[float]] = None,
        heights: Optional[Sequence[float]] = None,
        labels: Optional[Sequence[str]] = None,
    ):
        self.plots = list(plots)
        n = len(self.plots)
        if ncol is None and nrow is None:
            self.ncol = n if n <= 2 else int(math.ceil(math.sqrt(n)))
            self.nrow = int(math.ceil(n / max(self.ncol, 1)))
        elif ncol is not None:
            self.ncol = int(ncol)
            self.nrow = int(nrow) if nrow is not None else int(math.ceil(n / max(self.ncol, 1)))
        else:
            self.nrow = int(nrow)  # type: ignore[arg-type]
            self.ncol = int(math.ceil(n / max(self.nrow, 1)))
        self.widths = list(widths) if widths is not None else None
        self.heights = list(heights) if heights is not None else None
        self.labels = list(labels) if labels is not None else None

    def __len__(self) -> int:
        return len(self.plots)

    def __getitem__(self, idx: int) -> Any:
        return self.plots[idx]

    def __add__(self, other: Any) -> "PlotList":
        if isinstance(other, PlotModifier):
            # Apply modifier to the last plot in the list (just like `p + p + layout_dendrogram()` in R!)
            new_plots = list(self.plots)
            if new_plots:
                new_plots[-1] = new_plots[-1] + other
            return PlotList(new_plots, ncol=self.ncol, nrow=self.nrow, widths=self.widths, heights=self.heights, labels=self.labels)
        if isinstance(other, PlotList):
            return PlotList(self.plots + other.plots, ncol=self.ncol)
        return PlotList(self.plots + [other], ncol=self.ncol)

    def __truediv__(self, other: Any) -> "PlotList":
        if self.ncol == 1:
            extra = other.plots if isinstance(other, PlotList) and other.ncol == 1 else [other]
            all_p = self.plots + extra
            return PlotList(all_p, ncol=1, nrow=len(all_p))
        return PlotList([self, other], ncol=1, nrow=2)

    def __or__(self, other: Any) -> "PlotList":
        if self.nrow == 1:
            extra = other.plots if isinstance(other, PlotList) and other.nrow == 1 else [other]
            all_p = self.plots + extra
            return PlotList(all_p, ncol=len(all_p), nrow=1)
        return PlotList([self, other], ncol=2, nrow=1)

    def to_svg(self, width: int = 800, height: int = 600, _in_grid: bool = False) -> str:
        from .plot import render_plotlist_svg
        return render_plotlist_svg(self, width=width, height=height)

    def _render_rgba(self, width: int = 800, height: int = 600) -> np.ndarray:
        from .plot import render_plotlist_rgba
        return render_plotlist_rgba(self, width=width, height=height)

    def to_png(self, width: int = 800, height: int = 600) -> bytes:
        from .plot import render_plotlist_png
        return render_plotlist_png(self, width=width, height=height)

    def save(self, filename: str, width: float = 8.0, height: float = 6.0, dpi: int = 100) -> str:
        from .plot import ggsave
        return ggsave(filename, plot=self, width=width, height=height, dpi=dpi)

    def _repr_svg_(self) -> str:
        return self.to_svg()


def multiplot(*plots: Any, plotlist: Optional[Sequence[Any]] = None, ncol: int = 1, widths: Optional[Sequence[float]] = None, labels: Optional[Sequence[str]] = None) -> PlotList:
    """Place multiple plots on one page (exact port of ggtree::multiplot in R/multiplot.R)."""
    all_plots = list(plots) + (list(plotlist) if plotlist is not None else [])
    return PlotList(all_plots, ncol=ncol, widths=widths, labels=labels)


def plot_list(*plots: Any, ncol: Optional[int] = None, nrow: Optional[int] = None, widths: Optional[Sequence[float]] = None, heights: Optional[Sequence[float]] = None, tag_levels: Optional[str] = None) -> PlotList:
    """Combine multiple plots into a composite PlotList (compatible with aplot::plot_list)."""
    flat = []
    for p in plots:
        if isinstance(p, (list, tuple)) and not isinstance(p, PlotList):
            flat.extend(p)
        else:
            flat.append(p)
    labels = None
    if tag_levels == "A":
        labels = [chr(ord("A") + i) for i in range(len(flat))]
    elif tag_levels == "a":
        labels = [chr(ord("a") + i) for i in range(len(flat))]
    return PlotList(flat, ncol=ncol, nrow=nrow, widths=widths, heights=heights, labels=labels)
