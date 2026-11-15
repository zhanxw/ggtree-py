"""
Core GGTree / ggplot object, grammar-of-graphics composition, ggplot_build, and native SVG/PNG/PDF renderer.
"""

from __future__ import annotations
import copy
import html
import math
import re
import shutil
import struct
import subprocess
import zlib
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .dataframe import DataFrame, TblTree, NA, is_na, _is_scalar_na, getOption, ensure_dataframe, is_dataframe_like
from .tree import Phylo, TreeData, MultiPhylo, TreeDataList
from .layout import fortify, get_clade_position_, check_reverse
from .aes import Aes, aes, modifyList, eval_expr, eval_subset
from .themes_scales import (
    PlotModifier, Theme, theme_tree, Scale, Coord, coord_cartesian,
    coord_polar, coord_flip, revts
)
from .geoms import Layer, LayerList, geom_tree
from .export import (
    InteractivePlot,
    append_interactive_node_overlay,
    format_tile_bar_tooltip,
    render_interactive_html,
    svg_to_pdf_bytes,
)


_LAST_PLOT: Optional["GGTree"] = None


def set_last_plot(plot: "GGTree") -> None:
    global _LAST_PLOT
    _LAST_PLOT = plot


def last_plot() -> Optional["GGTree"]:
    return _LAST_PLOT


class _RangeBox:
    def __init__(self, r: Sequence[float]):
        self.range = list(r)

    def __getitem__(self, key: str) -> Any:
        if key == "range":
            return self.range
        raise KeyError(key)


class _ScaleBox:
    def __init__(self, r: Sequence[float]):
        self.range = _RangeBox(r)

    def __getitem__(self, key: str) -> Any:
        if key == "range":
            return self.range
        raise KeyError(key)


class BuiltLayout:
    """Simulated ggplot_build()$layout structure compatible with R ggtree tests."""

    def __init__(
        self,
        x_range: Tuple[float, float],
        y_range: Tuple[float, float],
        panels: List[str],
        panel_limits: Optional[Dict[str, Dict[str, Tuple[Optional[float], Optional[float]]]]] = None,
    ):
        panel_limits = panel_limits or {}
        self.panel_scales_x: List[_ScaleBox] = []
        self.panel_scales_y: List[_ScaleBox] = []
        self.panel_params: List[Dict[str, Any]] = []
        self.panels = panels

        for pan in (panels or ["Tree"]):
            px = list(x_range)
            py = list(y_range)
            if pan in panel_limits:
                if "x" in panel_limits[pan]:
                    lo, hi = panel_limits[pan]["x"]
                    if lo is not None:
                        px[0] = float(lo)
                    if hi is not None:
                        px[1] = float(hi)
                if "y" in panel_limits[pan]:
                    lo, hi = panel_limits[pan]["y"]
                    if lo is not None:
                        py[0] = float(lo)
                    if hi is not None:
                        py[1] = float(hi)
            self.panel_scales_x.append(_ScaleBox(px))
            self.panel_scales_y.append(_ScaleBox(py))
            self.panel_params.append({
                "x.range": px,
                "y.range": py,
                "x": {"range": px},
                "y": {"range": py},
            })


class BuiltPlot:
    """Result of ggplot_build(p), providing .data, .layout, and .plot."""

    def __init__(
        self,
        plot: "GGTree",
        layer_data: List[DataFrame],
        x_range: Tuple[float, float],
        y_range: Tuple[float, float],
        panels: List[str],
        panel_limits: Optional[Dict[str, Dict[str, Tuple[Optional[float], Optional[float]]]]] = None,
    ):
        self.plot = plot
        self.data = layer_data
        self.layout = BuiltLayout(x_range, y_range, panels, panel_limits=panel_limits)

    def __getitem__(self, key: str) -> Any:
        if key == "data":
            return self.data
        if key == "layout":
            return self.layout
        if key == "plot":
            return self.plot
        raise KeyError(key)


class GGTree:
    """
    Grammar-of-graphics phylogenetic tree plot object (compatible with R ggplot / ggtree).
    Supports:
    - `p + geom_*() + theme_*() + scale_*()`
    - `p %<+% df` via `p % df`, `p << df`, `p @ df`, or `p + df`
    - `p %<% tr2` via `p % tr2`, `p << tr2`, `p @ tr2`, or `p + tr2`
    - `p1 | p2` (side-by-side) and `p1 / p2` (stacked)
    - `ggplot_build(p)`
    - `ggsave(filename, p)` / `p.save(filename)` / `p.to_svg()` / `p.to_png()`
    """

    is_ggtree: bool = True

    def __init__(
        self,
        data: Optional[DataFrame] = None,
        mapping: Optional[Aes] = None,
        layout: str = "rectangular",
        tr: Any = None,
        extra_params: Optional[Dict[str, Any]] = None,
    ):
        self.data: DataFrame = data if data is not None else DataFrame()
        self.mapping: Aes = mapping.copy() if isinstance(mapping, Aes) else Aes(**(mapping or {}))
        self.layers: List[Layer] = []
        self.scales: Dict[str, Scale] = {}
        self.theme_spec: Dict[str, Any] = {"name": "tree", "bgcolor": "white", "fgcolor": "black", "show_x_axis": False}
        self.guides: Dict[str, Any] = {}
        self.coord: Coord = coord_cartesian()
        self.labels: Dict[str, Any] = {}
        self.layout: str = layout if isinstance(layout, str) else "rectangular"
        self._source_tr: Any = tr
        self._extra_params: Dict[str, Any] = extra_params or {}
        self.plot_env: Dict[str, Any] = {}
        self.facet_spec: Any = None
        self.panels: List[str] = ["Tree"]
        self.panel_labels: Dict[str, str] = {}
        self.panel_widths: Dict[Any, float] = {}
        self.panel_limits: Dict[str, Dict[str, Tuple[Optional[float], Optional[float]]]] = {}
        self.facet_data_store: Dict[str, DataFrame] = {}
        self.xlim: Optional[Tuple[Optional[float], Optional[float]]] = None
        self.ylim: Optional[Tuple[Optional[float], Optional[float]]] = None
        self.open_angle: float = 0.0
        self.rotate_angle: float = 0.0
        self.attrs: Dict[str, Any] = {}
        set_last_plot(self)

    @property
    def layout_name(self) -> str:
        return self.layout

    @layout_name.setter
    def layout_name(self, val: str) -> None:
        self.layout = val

    @property
    def _xlim(self) -> Optional[Tuple[Optional[float], Optional[float]]]:
        if self.xlim is not None:
            return self.xlim
        if "Tree" in self.panel_limits and "x" in self.panel_limits["Tree"]:
            return self.panel_limits["Tree"]["x"]
        return None

    @property
    def _ylim(self) -> Optional[Tuple[Optional[float], Optional[float]]]:
        if self.ylim is not None:
            return self.ylim
        if "Tree" in self.panel_limits and "y" in self.panel_limits["Tree"]:
            return self.panel_limits["Tree"]["y"]
        return None

    @property
    def panel_xlims(self) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
        res: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
        for pan, lims in self.panel_limits.items():
            if "x" in lims:
                res[pan] = lims["x"]
        if self.xlim is not None and "Tree" not in res:
            res["Tree"] = self.xlim
        return res

    def copy(self) -> "GGTree":
        new = GGTree(
            data=self.data.copy() if self.data is not None else DataFrame(),
            mapping=self.mapping.copy() if self.mapping is not None else Aes(),
            layout=self.layout,
            tr=self._source_tr,
            extra_params=dict(self._extra_params),
        )
        new.layers = [l.copy() for l in self.layers]
        new.scales = {k: copy.copy(v) for k, v in self.scales.items()}
        new.theme_spec = dict(self.theme_spec)
        new.guides = dict(self.guides)
        new.coord = copy.copy(self.coord)
        new.labels = dict(self.labels)
        new.plot_env = dict(self.plot_env)
        new.facet_spec = copy.copy(self.facet_spec)
        new.panels = list(self.panels)
        new.panel_labels = dict(self.panel_labels)
        new.panel_widths = dict(self.panel_widths)
        new.panel_limits = {k: dict(v) for k, v in self.panel_limits.items()}
        new.facet_data_store = {k: v.copy() for k, v in self.facet_data_store.items()}
        new.xlim = self.xlim
        new.ylim = self.ylim
        new.open_angle = self.open_angle
        new.rotate_angle = self.rotate_angle
        new.attrs = copy.deepcopy(self.attrs)
        set_last_plot(new)
        return new

    def __getitem__(self, key: str) -> Any:
        if key == "data":
            return self.data
        if key == "layers":
            return self.layers
        if key == "mapping":
            return self.mapping
        if key == "scales":
            return self.scales
        if key == "theme":
            return self.theme_spec
        if key in ("coordinates", "coord"):
            return self.coord
        if key == "labels":
            return self.labels
        if key == "plot_env":
            return self.plot_env
        raise KeyError(key)

    def __setitem__(self, key: str, value: Any) -> None:
        if key == "data":
            self.data = value
        elif key == "layers":
            self.layers = value
        elif key == "mapping":
            self.mapping = value
        elif key == "scales":
            self.scales = value
        elif key == "theme":
            self.theme_spec = value
        elif key in ("coordinates", "coord"):
            self.coord = value
        elif key == "labels":
            self.labels = value
        else:
            setattr(self, key, value)

    # ------------------------------------------------------------------
    # Data attachment (%<+%) and tree replacement (%<%)
    # ------------------------------------------------------------------
    def attach_data(self, df: Any) -> "GGTree":
        """Implements R's `p %<+% df`: left-joins taxon/node metadata into plot data."""
        df = ensure_dataframe(df)
        p = self.copy()
        if len(df) == 0 or len(df.columns) == 0:
            return p

        first_col = df.columns[0]
        first_vals = df[first_col].tolist()
        labels_set = set(str(v) for v in p.data["label"].tolist() if not _is_scalar_na(v)) if "label" in p.data else set()
        nodes_set = set(p.data["node"].tolist()) if "node" in p.data else set()

        match_label = sum(1 for v in first_vals if not _is_scalar_na(v) and str(v) in labels_set)
        match_node = sum(1 for v in first_vals if not _is_scalar_na(v) and v in nodes_set)

        join_source_col: str = first_col
        if "node" in df.columns and match_node > 0 and first_col == "node":
            key_col = "node"
        elif match_label >= match_node and match_label > 0:
            key_col = "label"
        elif match_node > 0:
            key_col = "node"
        elif "label" in df.columns and any(not _is_scalar_na(v) and str(v) in labels_set for v in df["label"].tolist()):
            key_col = "label"
            join_source_col = "label"
        elif "node" in df.columns and any(not _is_scalar_na(v) and v in nodes_set for v in df["node"].tolist()):
            key_col = "node"
            join_source_col = "node"
        elif df._rownames is not None and len(df._rownames) == len(df):
            rn_label_hits = sum(1 for rn in df._rownames if rn in labels_set)
            rn_node_vals = [int(rn) if rn.lstrip("-").isdigit() else rn for rn in df._rownames]
            rn_node_hits = sum(1 for rv in rn_node_vals if rv in nodes_set)
            if rn_label_hits > 0 and rn_label_hits >= rn_node_hits:
                key_col = "label"
                df = df.copy()
                df[".__rowname_key__"] = list(df._rownames)
                join_source_col = ".__rowname_key__"
            elif rn_node_hits > 0:
                key_col = "node"
                df = df.copy()
                df[".__rowname_key__"] = rn_node_vals
                join_source_col = ".__rowname_key__"
            else:
                key_col = "label" if match_label >= match_node else "node"
        else:
            key_col = "node"

        lookup: Dict[Any, Dict[str, Any]] = {}
        for r in df.to_records():
            k = r.get(join_source_col)
            if not _is_scalar_na(k):
                lookup[k] = r
                lookup[str(k)] = r

        p_records = p.data.to_records()
        for col in df.columns:
            if col in (key_col, join_source_col, ".__rowname_key__"):
                continue
            new_vals = []
            for r in p_records:
                kv = r.get(key_col)
                matched = lookup.get(kv)
                if matched is None and not _is_scalar_na(kv):
                    matched = lookup.get(str(kv))
                if matched is not None and col in matched:
                    new_vals.append(matched[col])
                else:
                    new_vals.append(r.get(col, None))
            p.data[col] = new_vals
            if col in df._levels:
                p.data._levels[col] = list(df._levels[col])

        set_last_plot(p)
        return p

    def replace_tree(self, tr: Any) -> "GGTree":
        """Implements R's `p %<% tr2`: replaces tree data while preserving all layers/styling."""
        p = self.copy()
        params = dict(self._extra_params)
        p._source_tr = tr
        p.data = fortify(
            tr,
            layout=p.layout,
            ladderize=params.get("ladderize", True),
            right=params.get("right", False),
            branch_length=params.get("branch_length", "branch.length"),
            mrsd=params.get("mrsd", None),
            as_Date=params.get("as_Date", False),
            yscale=params.get("yscale", "none"),
            yscale_mapping=params.get("yscale_mapping", None),
            root_position=params.get("root_position", 0.0),
        )
        set_last_plot(p)
        return p

    # ------------------------------------------------------------------
    # Operators: +, %, <<, @, |, /
    # ------------------------------------------------------------------
    def __add__(self, other: Any) -> Any:
        if other is None:
            return self
        if isinstance(other, (list, LayerList)):
            p = self.copy()
            for item in other:
                p = p + item
            set_last_plot(p)
            return p
        if isinstance(other, Layer):
            p = other.apply(self)
            set_last_plot(p)
            return p
        if isinstance(other, PlotModifier):
            p = other.apply(self.copy())
            set_last_plot(p)
            return p
        if isinstance(other, Aes):
            p = self.copy()
            p.mapping = modifyList(p.mapping, other)
            set_last_plot(p)
            return p
        if isinstance(other, (Phylo, TreeData, MultiPhylo, TreeDataList)):
            return self.replace_tree(other)
        if is_dataframe_like(other):
            return self.attach_data(other)
        if isinstance(other, GGTree):
            from .compose import PlotList
            return PlotList([self.copy(), other.copy()])
        if callable(other):
            res = other(self.copy())
            set_last_plot(res)
            return res
        raise TypeError(f"Cannot add object of type {type(other).__name__} to GGTree")

    def __mod__(self, other: Any) -> "GGTree":
        """Supports `p % df` (like `p %<+% df`) and `p % tr2` (like `p %<% tr2`)."""
        if isinstance(other, (Phylo, TreeData, MultiPhylo, TreeDataList)):
            return self.replace_tree(other)
        if isinstance(other, GGTree):
            p = self.copy()
            for layer in other.layers:
                if layer.geom_type not in ("tree", "tree_horizontal", "tree_vertical"):
                    p = p + layer
            return p
        return self.attach_data(other)

    def __lshift__(self, other: Any) -> "GGTree":
        return self.__mod__(other)

    def __matmul__(self, other: Any) -> "GGTree":
        return self.__mod__(other)

    def __or__(self, other: Any) -> Any:
        """Horizontal side-by-side composition (`p1 | p2`)."""
        from .compose import PlotList
        if isinstance(other, PlotList):
            return PlotList([self] + other.plots, ncol=len(other.plots) + 1)
        if isinstance(other, GGTree):
            return PlotList([self, other], ncol=2, nrow=1)
        raise TypeError(f"Unsupported operand for |: {type(other).__name__}")

    def __truediv__(self, other: Any) -> Any:
        """Vertical stacked composition (`p1 / p2`)."""
        from .compose import PlotList
        if isinstance(other, PlotList):
            return PlotList([self] + other.plots, ncol=1, nrow=len(other.plots) + 1)
        if isinstance(other, GGTree):
            return PlotList([self, other], ncol=1, nrow=2)
        raise TypeError(f"Unsupported operand for /: {type(other).__name__}")

    # ------------------------------------------------------------------
    # Building plot layers & ranges (ggplot_build)
    # ------------------------------------------------------------------
    def build(self) -> BuiltPlot:
        layer_dfs = [layer.compute_data(self.data, self.mapping, layout=self.layout) for layer in self.layers]

        xs: List[float] = []
        ys: List[float] = []
        for df in [self.data] + layer_dfs:
            for col in ("x", "xmin", "xmax", "xend", "branch"):
                if col in df.columns:
                    for v in df[col].tolist():
                        if isinstance(v, (int, float, np.number)) and not _is_scalar_na(v):
                            xs.append(float(v))
            for col in ("y", "ymin", "ymax", "yend"):
                if col in df.columns:
                    for v in df[col].tolist():
                        if isinstance(v, (int, float, np.number)) and not _is_scalar_na(v):
                            ys.append(float(v))

        xmin = min(xs) if xs else 0.0
        xmax = max(xs) if xs else 1.0
        ymin = min(ys) if ys else 0.0
        ymax = max(ys) if ys else 1.0

        if self.xlim is not None:
            if self.xlim[0] is not None:
                xmin = float(self.xlim[0])
            if self.xlim[1] is not None:
                xmax = float(self.xlim[1])
        if self.ylim is not None:
            if self.ylim[0] is not None:
                ymin = float(self.ylim[0])
            if self.ylim[1] is not None:
                ymax = float(self.ylim[1])

        exp_x = (xmax - xmin) * 0.05 if xmax > xmin else 0.05
        exp_y = (ymax - ymin) * 0.05 if ymax > ymin else 0.05
        x_range = (xmin if (self.xlim and self.xlim[0] is not None) else xmin - exp_x,
                   xmax if (self.xlim and self.xlim[1] is not None) else xmax + exp_x)
        y_range = (ymin if (self.ylim and self.ylim[0] is not None) else ymin - exp_y,
                   ymax if (self.ylim and self.ylim[1] is not None) else ymax + exp_y)

        panels = list(self.panels) if self.panels else ["Tree"]
        return BuiltPlot(self, layer_dfs, x_range, y_range, panels, panel_limits=self.panel_limits)

    # ------------------------------------------------------------------
    # Native SVG / PNG rendering (zero external dependencies)
    # ------------------------------------------------------------------
    def to_svg(
        self,
        width: float = 800,
        height: float = 600,
        _in_grid: bool = False,
        _interactive_overlay: bool = False,
    ) -> str:
        """Render the GGTree plot to a standalone SVG string."""
        built = self.build()
        w = max(int(width), 100)
        h = max(int(height), 80)

        x_min, x_max = built.layout.panel_scales_x[0].range.range
        y_min, y_max = built.layout.panel_scales_y[0].range.range
        dx = (x_max - x_min) if x_max != x_min else 1.0
        dy = (y_max - y_min) if y_max != y_min else 1.0

        is_polar = self.coord.coord_type == "polar" or self.layout in ("circular", "fan", "inward_circular", "radial")
        is_flipped = self.coord.coord_type == "flip" or self.layout == "dendrogram"
        is_multi_panel = len(self.panels) > 1 and not is_polar

        # Inspect layers for dynamic margins, axis, legends, and text extent
        show_x_axis = bool(self.theme_spec.get("show_x_axis", False))
        show_y_axis = bool(self.theme_spec.get("show_y_axis", False))
        has_polar_labels = False
        has_rotated_bottom_text = False
        max_right_label_chars = 0
        has_tile_legend = False
        has_color_scale_legend = False

        for layer, ldf in zip(self.layers, built.data):
            gt = layer.geom_type
            if gt == "tile" and len(ldf) > 0:
                has_tile_legend = True
            inherit_flag = layer.params.get("inherit.aes", layer.params.get("inherit_aes", True))
            merged_aes = self.mapping.merge(layer.mapping) if inherit_flag else layer.mapping
            if gt not in ("tile", "blank", "inset") and len(ldf) > 0:
                if "color" in merged_aes or "colour" in merged_aes or "fill" in merged_aes:
                    c_key = "fill" if "fill" in merged_aes else ("color" if "color" in merged_aes else "colour")
                    scale_hint = layer.params.get(f"_{c_key}_scale") or self.scales.get("fill" if c_key == "fill" else "color")
                    col_vals = [r.get(c_key, r.get("color", r.get("colour", r.get("fill")))) for r in ldf.to_records()]
                    non_na_vals = [v for v in col_vals if not _is_scalar_na(v) and str(v) != ""]
                    if non_na_vals and (scale_hint is not None or not all(_is_direct_css_color(v) for v in non_na_vals)):
                        has_color_scale_legend = True
            if gt in ("tiplab", "tiplab2", "nodelab", "nodelab2", "text", "text2", "label", "label2", "cladelab_text", "strip_text") and len(ldf) > 0:
                ang_param = layer.params.get("angle", 0.0)
                if isinstance(ang_param, (int, float)) and abs(float(ang_param)) > 15:
                    has_rotated_bottom_text = True
                for r in ldf.to_records():
                    lbl = r.get("label")
                    if _is_scalar_na(lbl) or lbl is None:
                        continue
                    s_lbl, _, _ = _parse_plotmath_label(str(lbl), bool(r.get("parse", layer.params.get("parse", False))))
                    if not s_lbl:
                        continue
                    if is_polar:
                        has_polar_labels = True
                    max_right_label_chars = max(max_right_label_chars, len(s_lbl))

        compact = _in_grid or w < 460 or h < 340
        has_facet_strips = is_multi_panel
        title = self.labels.get("title")

        pad_top = (36 if title else 18) if not compact else (24 if title else 10)
        if has_facet_strips:
            pad_top += 24
        pad_bottom = 20 if compact else 30
        if show_x_axis or is_multi_panel:
            pad_bottom = max(pad_bottom, 42 if not compact else 28)
        if has_rotated_bottom_text and not compact:
            pad_bottom = max(pad_bottom, 74)
        if is_flipped and max_right_label_chars > 0 and not compact:
            pad_bottom = max(pad_bottom, 58)

        pad_left = (38 if show_y_axis else 18) if compact else (52 if show_y_axis else 28)
        pad_right = 18 if compact else 32
        if not is_polar and not is_flipped and max_right_label_chars > 0 and not compact:
            # Reserve extra space when non-facet right-hand labels sit near x_max
            est_px = min(145, max(36, int(max_right_label_chars * 6.2)))
            pad_right = max(pad_right, est_px)
        if has_tile_legend and not compact:
            pad_right = max(pad_right, 118)
        elif has_color_scale_legend and not compact:
            if not is_polar and not is_multi_panel and max_right_label_chars > 0:
                pad_right = max(pad_right + 88, 128)
            else:
                pad_right = max(pad_right, 118)

        plot_w = max(w - pad_left - pad_right, 60)
        plot_h = max(h - pad_top - pad_bottom, 60)

        # Polar coordinate calibration using tip y span so tip 1 and tip N don't overlap at 0/360
        tip_ys = [
            float(r["y"])
            for r in self.data.to_records()
            if r.get("isTip") and not _is_scalar_na(r.get("y"))
        ]
        if not tip_ys:
            tip_ys = [float(v) for v in self.data["y"].tolist() if not _is_scalar_na(v)] if "y" in self.data else [y_min, y_max]
        tip_y_lo = min(tip_ys) if tip_ys else y_min
        tip_y_hi = max(tip_ys) if tip_ys else y_max
        open_ang = float(getattr(self, "open_angle", 0.0) or 0.0)
        rot_ang = float(getattr(self, "rotate_angle", 0.0) or 0.0)
        if self.layout == "fan" and abs(open_ang) < 1e-6:
            open_ang = 45.0

        def polar_deg(y: float) -> float:
            if open_ang > 1e-6:
                span_deg = max(360.0 - open_ang, 10.0)
                t = (y - tip_y_lo) / max(tip_y_hi - tip_y_lo, 1.0)
                return t * span_deg + rot_ang
            t = (y - tip_y_lo + 0.5) / max((tip_y_hi - tip_y_lo) + 1.0, 1.0)
            return t * 360.0 + rot_ang

        cx, cy = pad_left + plot_w / 2.0, pad_top + plot_h / 2.0
        max_r_factor = 0.31 if has_polar_labels else 0.43
        if compact:
            max_r_factor = 0.36 if has_polar_labels else 0.43
        max_r = min(plot_w, plot_h) * max_r_factor
        inner_r = max_r * (0.24 if self.layout == "inward_circular" else (0.04 if self.layout == "circular" else 0.0))

        def polar_rad(x: float) -> float:
            r_norm = (x - x_min) / dx if dx != 0 else 0.0
            if self.layout == "inward_circular":
                return max(0.0, max_r - r_norm * (max_r - inner_r))
            return max(0.0, inner_r + r_norm * (max_r - inner_r))

        # Panel geometry computation (single panel vs facet_plot multi-panel)
        panel_boxes: Dict[str, Tuple[float, float, float, float, float, float]] = {}
        if is_multi_panel:
            gap = 16.0
            n_pan = len(self.panels)
            raw_ws = [float(self.panel_widths.get(pan, 1.0)) for pan in self.panels]
            sum_ws = sum(max(v, 0.05) for v in raw_ws) or 1.0
            avail_w = max(plot_w - gap * (n_pan - 1), 40.0)
            cur_x = float(pad_left)
            for pan, rw in zip(self.panels, raw_ws):
                pw = avail_w * (max(rw, 0.05) / sum_ws)
                # Compute per-panel x limits from non-blank layers, overridden by panel_limits
                pan_xs: List[float] = []
                if pan == "Tree" and "x" in self.data:
                    pan_xs.extend(float(v) for v in self.data["x"].tolist() if not _is_scalar_na(v))
                for layer, ldf in zip(self.layers, built.data):
                    if layer.geom_type == "blank":
                        continue
                    layer_pan = layer.params.get("panel")
                    has_pcol = "PANEL" in ldf.columns and len(ldf) > 0
                    in_pan = (layer_pan == pan) or (not layer_pan and (not has_pcol or str(ldf["PANEL"][0]) == pan))
                    if in_pan:
                        for col in ("x", "xmin", "xmax", "xend", "branch"):
                            if col in ldf.columns:
                                for v in ldf[col].tolist():
                                    if isinstance(v, (int, float, np.number)) and not _is_scalar_na(v):
                                        pan_xs.append(float(v))
                        if layer.geom_type in ("col", "bar") and len(ldf) > 0:
                            pan_xs.append(0.0)
                px_lo = min(pan_xs) if pan_xs else x_min
                px_hi = max(pan_xs) if pan_xs else x_max
                lim_x = self.panel_limits.get(pan, {}).get("x")
                span_p = (px_hi - px_lo) if px_hi > px_lo else 1.0
                px0 = float(lim_x[0]) if (lim_x and lim_x[0] is not None) else (px_lo - 0.05 * span_p)
                px1 = float(lim_x[1]) if (lim_x and lim_x[1] is not None) else (px_hi + 0.05 * span_p)
                panel_boxes[pan] = (cur_x, float(pad_top), pw, float(plot_h), px0, px1)
                cur_x += pw + gap
        else:
            panel_boxes["Tree"] = (float(pad_left), float(pad_top), float(plot_w), float(plot_h), float(x_min), float(x_max))

        def map_pt(x: float, y: float, panel: str = "Tree") -> Tuple[float, float]:
            if is_polar:
                r = polar_rad(x)
                theta = math.radians(polar_deg(y))
                return (cx + r * math.cos(theta), cy - r * math.sin(theta))
            bx, by, bw, bh, px0, px1 = panel_boxes.get(panel, panel_boxes["Tree"])
            pdx = (px1 - px0) if px1 != px0 else 1.0
            if is_flipped:
                px = bx + ((y - y_min) / dy) * bw
                py = by + ((x - px0) / pdx) * bh
                return (px, py)
            px = bx + ((x - px0) / pdx) * bw
            py = by + bh - ((y - y_min) / dy) * bh
            return (px, py)

        def polar_arc_d(x_const: float, y0: float, y1: float, move_cmd: str = "M") -> str:
            r = polar_rad(x_const)
            d0 = polar_deg(y0)
            d1 = polar_deg(y1)
            n_steps = max(2, int(math.ceil(abs(d1 - d0) / 3.0)))
            seg_pts: List[str] = []
            for k in range(n_steps + 1):
                dk = d0 + (d1 - d0) * (k / n_steps)
                th = math.radians(dk)
                px = cx + r * math.cos(th)
                py = cy - r * math.sin(th)
                cmd = move_cmd if k == 0 else "L"
                seg_pts.append(f"{cmd} {px:.2f} {py:.2f}")
            return " ".join(seg_pts)

        def polar_wedge_d(r_x0: float, r_x1: float, y0: float, y1: float) -> str:
            r_in = polar_rad(r_x0)
            r_out = polar_rad(r_x1)
            d0 = polar_deg(y0)
            d1 = polar_deg(y1)
            n_steps = max(2, int(math.ceil(abs(d1 - d0) / 3.0)))
            cmds: List[str] = []
            for k in range(n_steps + 1):
                th = math.radians(d0 + (d1 - d0) * (k / n_steps))
                px = cx + r_out * math.cos(th)
                py = cy - r_out * math.sin(th)
                cmds.append(f"{'M' if k == 0 else 'L'} {px:.2f} {py:.2f}")
            for k in range(n_steps, -1, -1):
                th = math.radians(d0 + (d1 - d0) * (k / n_steps))
                px = cx + r_in * math.cos(th)
                py = cy - r_in * math.sin(th)
                cmds.append(f"L {px:.2f} {py:.2f}")
            cmds.append("Z")
            return " ".join(cmds)

        bg = str(self.theme_spec.get("bgcolor", "white") or "white")
        if bg == "transparent":
            bg = "white"
        parts: List[str] = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
            f'<rect width="{w}" height="{h}" fill="{html.escape(bg, quote=True)}"/>',
        ]

        if title:
            t_fsize = 12 if compact else 15
            t_y = 17 if compact else 24
            parts.append(
                f'<text x="{w / 2:.1f}" y="{t_y}" text-anchor="middle" font-family="sans-serif" font-size="{t_fsize}" font-weight="bold" fill="#1F2937">{html.escape(str(title))}</text>'
            )

        # Panel backgrounds, facet header strips, and subtle gridlines
        if is_multi_panel:
            for pan in self.panels:
                bx, by, bw, bh, px0, px1 = panel_boxes[pan]
                pan_lbl = str(self.panel_labels.get(pan, pan))
                strip_h = 20.0
                parts.append(
                    f'<rect x="{bx:.2f}" y="{by - strip_h:.2f}" width="{bw:.2f}" height="{strip_h:.2f}" fill="#E5E7EB" stroke="#9CA3AF" stroke-width="0.8"/>'
                )
                parts.append(
                    f'<text x="{bx + bw / 2:.2f}" y="{by - 6:.2f}" text-anchor="middle" font-family="sans-serif" font-size="11" font-weight="bold" fill="#111827">{html.escape(pan_lbl)}</text>'
                )
                parts.append(
                    f'<rect x="{bx:.2f}" y="{by:.2f}" width="{bw:.2f}" height="{bh:.2f}" fill="#FAFAFA" stroke="#D1D5DB" stroke-width="0.8"/>'
                )
                # Bottom ticks per panel
                x_lbl_fn = self.scales["x"].params.get("labels") if "x" in self.scales else None
                for tv, tlbl in _nice_ticks(px0, px1, n=4, label_fn=x_lbl_fn):
                    tx, _ = map_pt(tv, y_min, panel=pan)
                    parts.append(f'<line x1="{tx:.2f}" y1="{by + bh:.2f}" x2="{tx:.2f}" y2="{by + bh + 4:.2f}" stroke="#4B5563" stroke-width="1"/>')
                    parts.append(f'<text x="{tx:.2f}" y="{by + bh + 15:.2f}" text-anchor="middle" font-family="sans-serif" font-size="9.5" fill="#4B5563">{html.escape(tlbl)}</text>')
        elif not is_polar and str(self.theme_spec.get("name", "tree") or "tree") in ("bw", "minimal"):
            bx, by, bw, bh, px0, px1 = panel_boxes["Tree"]
            for tv, _ in _nice_ticks(px0, px1, n=5):
                tx, _ = map_pt(tv, y_min)
                parts.append(f'<line x1="{tx:.2f}" y1="{by:.2f}" x2="{tx:.2f}" y2="{by + bh:.2f}" stroke="#E5E7EB" stroke-width="0.8"/>')
            for yv, _ in _nice_ticks(y_min, y_max, n=5):
                _, ty = map_pt(px0, yv)
                parts.append(f'<line x1="{bx:.2f}" y1="{ty:.2f}" x2="{bx + bw:.2f}" y2="{ty:.2f}" stroke="#E5E7EB" stroke-width="0.8"/>')
            if str(self.theme_spec.get("name", "")) == "bw":
                parts.append(f'<rect x="{bx:.2f}" y="{by:.2f}" width="{bw:.2f}" height="{bh:.2f}" fill="none" stroke="#374151" stroke-width="1.0"/>')

        legends_to_draw: List[Dict[str, Any]] = []

        # Helper to resolve parent->child bezier curves cleanly for roundrect / ellipse
        if self.layout in ("roundrect", "ellipse"):
            recs_t = self.data.to_records()
            nmap_t = {r.get("node"): r for r in recs_t if r.get("node") is not None}
            tree_color = "#222222"
            tree_lw = 1.2 if not compact else 1.0
            for layer, ldf in zip(self.layers, built.data):
                if layer.geom_type in ("tree", "tree_horizontal"):
                    c0 = layer.params.get("color", layer.params.get("colour"))
                    if c0:
                        tree_color = str(c0)
                    lw0 = layer.params.get("linewidth", layer.params.get("size"))
                    if lw0 is not None:
                        tree_lw = max(float(lw0) * 1.4, 1.0)
                    break
            for r in recs_t:
                nid, pid = r.get("node"), r.get("parent")
                if nid is None or pid is None or nid == pid or pid not in nmap_t:
                    continue
                pr = nmap_t[pid]
                if any(_is_scalar_na(v) for v in (r.get("x"), r.get("y"), pr.get("x"), pr.get("y"))):
                    continue
                p_par = map_pt(float(pr["x"]), float(pr["y"]))
                p_cor = map_pt(float(pr["x"]), float(r["y"]))
                p_chi = map_pt(float(r["x"]), float(r["y"]))
                if self.layout == "ellipse":
                    parts.append(
                        f'<path d="M {p_par[0]:.2f} {p_par[1]:.2f} Q {p_cor[0]:.2f} {p_cor[1]:.2f} {p_chi[0]:.2f} {p_chi[1]:.2f}" fill="none" stroke="{html.escape(tree_color, quote=True)}" stroke-width="{tree_lw:.2f}"/>'
                    )
                else:
                    # Rounded rectangular elbow
                    vy = p_cor[1] - p_par[1]
                    hx = p_chi[0] - p_cor[0]
                    rad = min(abs(vy) * 0.45, abs(hx) * 0.45, 7.0 if not compact else 4.0)
                    if rad > 0.8:
                        sy = 1.0 if vy >= 0 else -1.0
                        sx = 1.0 if hx >= 0 else -1.0
                        c1y = p_cor[1] - sy * rad
                        c2x = p_cor[0] + sx * rad
                        parts.append(
                            f'<path d="M {p_par[0]:.2f} {p_par[1]:.2f} L {p_cor[0]:.2f} {c1y:.2f} Q {p_cor[0]:.2f} {p_cor[1]:.2f} {c2x:.2f} {p_cor[1]:.2f} L {p_chi[0]:.2f} {p_chi[1]:.2f}" fill="none" stroke="{html.escape(tree_color, quote=True)}" stroke-width="{tree_lw:.2f}"/>'
                        )
                    else:
                        parts.append(
                            f'<polyline points="{p_par[0]:.2f},{p_par[1]:.2f} {p_cor[0]:.2f},{p_cor[1]:.2f} {p_chi[0]:.2f},{p_chi[1]:.2f}" fill="none" stroke="{html.escape(tree_color, quote=True)}" stroke-width="{tree_lw:.2f}"/>'
                        )

        # Render each layer in z-order
        for layer, ldf in zip(self.layers, built.data):
            gt = layer.geom_type
            if gt == "blank":
                continue
            if self.layout in ("roundrect", "ellipse") and gt in ("tree", "tree2", "tree_horizontal", "tree_vertical"):
                continue

            raw_col = layer.params.get("color", layer.params.get("colour"))
            if isinstance(raw_col, (list, tuple)):
                color = str(raw_col[0]) if raw_col else "#222222"
            else:
                color = str(raw_col if raw_col is not None else "#222222")
            raw_fill = layer.params.get("fill")
            if isinstance(raw_fill, (list, tuple)):
                fill = str(raw_fill[0]) if raw_fill else "#4682B4"
            else:
                fill = str(raw_fill if raw_fill is not None else "#4682B4")
            alpha = float(layer.params.get("alpha", 1.0) if layer.params.get("alpha") is not None else 1.0)
            size = float(layer.params.get("size", layer.params.get("linewidth", 1.1)) or 1.1)
            default_pan = str(layer.params.get("panel", "Tree") or "Tree")
            if default_pan not in panel_boxes:
                default_pan = "Tree"
            inherit_flag = layer.params.get("inherit.aes", layer.params.get("inherit_aes", True))
            merged_aes = self.mapping.merge(layer.mapping) if inherit_flag else layer.mapping

            # 1. Highlight / Balance / Rectangles
            if gt in ("hilight", "highlight", "balance", "rect") and len(ldf) > 0:
                for r in ldf.to_records():
                    pan = str(r.get("PANEL", default_pan) or default_pan)
                    if pan not in panel_boxes:
                        pan = default_pan
                    x0 = float(r.get("xmin", r.get("x", 0.0)) or 0.0)
                    x1 = float(r.get("xmax", r.get("x", 0.0)) or 0.0)
                    y0 = float(r.get("ymin", r.get("y", 0.0)) or 0.0)
                    y1 = float(r.get("ymax", r.get("y", 0.0)) or 0.0)
                    fcol = _resolve_color_value(r.get("fill"), fill)
                    falpha = float(r.get("alpha", alpha if alpha < 1.0 else 0.35) or 0.35)
                    stroke_c = _resolve_color_value(r.get("color"), "none")
                    if is_polar:
                        d_attr = polar_wedge_d(x0, x1, y0 - 0.35, y1 + 0.35)
                        parts.append(
                            f'<path d="{d_attr}" fill="{html.escape(fcol, quote=True)}" fill-opacity="{falpha:.2f}" stroke="{html.escape(stroke_c, quote=True)}" stroke-width="0.6"/>'
                        )
                    else:
                        p0 = map_pt(x0, y1 + 0.35, panel=pan)
                        p1 = map_pt(x1, y0 - 0.35, panel=pan)
                        rx, ry = min(p0[0], p1[0]), min(p0[1], p1[1])
                        rw, rh = max(abs(p1[0] - p0[0]), 2.0), max(abs(p1[1] - p0[1]), 4.0)
                        parts.append(
                            f'<rect x="{rx:.2f}" y="{ry:.2f}" width="{rw:.2f}" height="{rh:.2f}" fill="{html.escape(fcol, quote=True)}" fill-opacity="{falpha:.2f}" stroke="{html.escape(stroke_c, quote=True)}" stroke-width="0.8"/>'
                        )

            # 2. Heatmap / MSA Tiles
            elif gt == "tile" and len(ldf) > 0:
                recs = ldf.to_records()
                has_val_mapping = ("value" in ldf.columns) or ("base" in ldf.columns) or ("_fill_scale" in layer.params) or ("fill" in layer.mapping)
                if has_val_mapping:
                    vals = [r.get("value", r.get("base", r.get("fill"))) for r in recs]
                    scale_obj = layer.params.get("_fill_scale") or (None if self.scales.get("fill") is None else self.scales.get("fill"))
                    is_msa = "base" in ldf.columns or layer.mapping.get("fill") == "base"
                    color_fn, leg_info = _build_value_color_mapper(
                        vals,
                        scale_obj,
                        default_color=fill,
                        is_msa=is_msa,
                        levels=(ldf._levels.get("value") or ldf._levels.get("fill")),
                    )
                    _add_legend_unique(legends_to_draw, leg_info)
                else:
                    color_fn = lambda v: _resolve_color_value(v, fill)

                border_c = str(layer.params.get("color", "white") or "white")
                for r in recs:
                    pan = str(r.get("PANEL", default_pan) or default_pan)
                    if pan not in panel_boxes:
                        pan = default_pan
                    if _is_scalar_na(r.get("x")) or _is_scalar_na(r.get("y")):
                        continue
                    tw = float(r.get("width", layer.params.get("width", dx * 0.05)) or dx * 0.05)
                    th = float(r.get("height", layer.params.get("height", 1.0)) or 1.0)
                    xc, yc = float(r["x"]), float(r["y"])
                    v_raw = r.get("value", r.get("base", r.get("fill")))
                    fcol = color_fn(v_raw)
                    if is_polar:
                        d_attr = polar_wedge_d(xc - tw / 2.0, xc + tw / 2.0, yc - th / 2.0, yc + th / 2.0)
                        tt_tile = format_tile_bar_tooltip(r, kind="tile") if _interactive_overlay else ""
                        parts.append(
                            f'<path d="{d_attr}" fill="{html.escape(fcol, quote=True)}" stroke="{html.escape(border_c, quote=True)}" stroke-width="0.4"{tt_tile}/>'
                        )
                    else:
                        p0 = map_pt(xc - tw / 2.0, yc + th / 2.0, panel=pan)
                        p1 = map_pt(xc + tw / 2.0, yc - th / 2.0, panel=pan)
                        rx, ry = min(p0[0], p1[0]), min(p0[1], p1[1])
                        rw, rh = max(abs(p1[0] - p0[0]), 1.2), max(abs(p1[1] - p0[1]), 1.2)
                        sw_tile = "0.2" if len(recs) > 250 else "0.5"
                        tt_tile = format_tile_bar_tooltip(r, kind="tile") if _interactive_overlay else ""
                        parts.append(
                            f'<rect x="{rx:.2f}" y="{ry:.2f}" width="{rw:.2f}" height="{rh:.2f}" fill="{html.escape(fcol, quote=True)}" stroke="{html.escape(border_c, quote=True)}" stroke-width="{sw_tile}"{tt_tile}/>'
                        )

            # 2b. Horizontal/facet column & bar charts (geom_col / geom_bar)
            elif gt in ("col", "bar") and len(ldf) > 0:
                recs_b = ldf.to_records()
                has_f_map = "fill" in merged_aes
                fmap_fn = None
                if has_f_map:
                    fvals = [r.get("fill") for r in recs_b]
                    fmap_fn, leg_b = _build_value_color_mapper(
                        fvals,
                        layer.params.get("_fill_scale") or self.scales.get("fill"),
                        default_color=fill,
                        force_discrete=True,
                        levels=(ldf._levels.get("fill") or self.data._levels.get(str(merged_aes.get("fill", "")))),
                    )
                    _add_legend_unique(legends_to_draw, leg_b)
                border_b = str(layer.params.get("color", layer.params.get("colour", "white")) or "white")
                for r in recs_b:
                    pan = str(r.get("PANEL", default_pan) or default_pan)
                    if pan not in panel_boxes:
                        pan = default_pan
                    if _is_scalar_na(r.get("x")) or _is_scalar_na(r.get("y")):
                        continue
                    xv = float(r["x"])
                    yc = float(r["y"])
                    bh_u = float(r.get("width", layer.params.get("width", 0.68)) or 0.68)
                    x_lo, x_hi = (0.0, xv) if xv >= 0.0 else (xv, 0.0)
                    p0 = map_pt(x_lo, yc + bh_u / 2.0, panel=pan)
                    p1 = map_pt(x_hi, yc - bh_u / 2.0, panel=pan)
                    rx, ry = min(p0[0], p1[0]), min(p0[1], p1[1])
                    rw, rh = max(abs(p1[0] - p0[0]), 1.2), max(abs(p1[1] - p0[1]), 1.5)
                    fcol = fmap_fn(r.get("fill")) if fmap_fn else _resolve_color_value(r.get("fill"), fill)
                    tt_bar = format_tile_bar_tooltip(r, kind="col") if _interactive_overlay else ""
                    parts.append(
                        f'<rect x="{rx:.2f}" y="{ry:.2f}" width="{rw:.2f}" height="{rh:.2f}" fill="{html.escape(fcol, quote=True)}" fill-opacity="{alpha:.2f}" stroke="{html.escape(border_b, quote=True)}" stroke-width="0.4"{tt_bar}/>'
                    )

            # 3. Polygons (collapsed clade triangles, etc.)
            elif gt == "polygon" and len(ldf) > 0:
                pts = []
                poly_fill = fill
                for r in ldf.to_records():
                    if not _is_scalar_na(r.get("fill")):
                        poly_fill = _resolve_color_value(r.get("fill"), poly_fill)
                    if not (_is_scalar_na(r.get("x")) or _is_scalar_na(r.get("y"))):
                        pan = str(r.get("PANEL", default_pan) or default_pan)
                        px, py = map_pt(float(r["x"]), float(r["y"]), panel=pan if pan in panel_boxes else default_pan)
                        pts.append(f"{px:.2f},{py:.2f}")
                if pts:
                    p_alpha = alpha if alpha < 1.0 else 0.65
                    parts.append(
                        f'<polygon points="{" ".join(pts)}" fill="{html.escape(poly_fill, quote=True)}" fill-opacity="{p_alpha:.2f}" stroke="{html.escape(color, quote=True)}" stroke-width="{max(size, 1.0):.2f}"/>'
                    )

            # 4. Line Segments & Tree Branches
            elif gt in ("tree", "tree2", "tree_horizontal", "tree_vertical", "segment", "segment2", "aline", "rootedge", "treescale", "cladelab_bar", "strip_bar", "taxalink", "tanglelink", "range") and len(ldf) > 0:
                lt_default = str(layer.params.get("linetype", "") or "")
                recs4 = ldf.to_records()
                cmap4 = None
                if ("color" in merged_aes or "colour" in merged_aes) and recs4:
                    c_expr4 = merged_aes.get("color", merged_aes.get("colour"))
                    cvals4 = [r.get("color", r.get("colour")) for r in recs4]
                    cmap4, leg4 = _build_value_color_mapper(
                        cvals4,
                        layer.params.get("_color_scale") or self.scales.get("color"),
                        default_color=color,
                        force_discrete=(str(c_expr4) in ("group", "focus")),
                        levels=(ldf._levels.get("color") or ldf._levels.get("colour") or self.data._levels.get(str(c_expr4))),
                    )
                    _add_legend_unique(legends_to_draw, leg4)
                for r in recs4:
                    pan = str(r.get("PANEL", default_pan) or default_pan)
                    if pan not in panel_boxes:
                        pan = default_pan
                    if any(_is_scalar_na(r.get(k)) for k in ("x", "y", "xend", "yend")):
                        continue
                    x0, y0, x1, y1 = float(r["x"]), float(r["y"]), float(r["xend"]), float(r["yend"])
                    if not all(math.isfinite(v) for v in (x0, y0, x1, y1)):
                        continue
                    rc = cmap4(r.get("color", r.get("colour"))) if cmap4 else _resolve_color_value(r.get("color", r.get("colour")), color)
                    rs = float(r.get("linewidth", r.get("size", size)) or size)
                    sw = max(rs * (1.8 if gt == "range" else 1.4), 1.0)
                    r_alpha = float(r.get("alpha", alpha) if r.get("alpha") is not None else alpha)
                    lt = str(r.get("linetype", lt_default) or "")
                    dash_attr = ' stroke-dasharray="2,3"' if lt == "dotted" else (' stroke-dasharray="5,3"' if lt == "dashed" else "")

                    if is_polar and gt in ("tree_vertical", "tree", "cladelab_bar", "strip_bar") and abs(x0 - x1) <= 1e-7 * max(dx, 1.0) and abs(y0 - y1) > 1e-7:
                        d_arc = polar_arc_d(x0, y0, y1, move_cmd="M")
                        parts.append(
                            f'<path d="{d_arc}" fill="none" stroke="{html.escape(rc, quote=True)}" stroke-width="{sw:.2f}" stroke-opacity="{r_alpha:.2f}"{dash_attr}/>'
                        )
                        continue

                    p0 = map_pt(x0, y0, panel=pan)
                    p1 = map_pt(x1, y1, panel=pan)
                    if abs(p0[0] - p1[0]) < 1e-6 and abs(p0[1] - p1[1]) < 1e-6:
                        continue
                    if gt == "taxalink":
                        if is_polar:
                            parts.append(
                                f'<path d="M {p0[0]:.2f} {p0[1]:.2f} Q {cx:.2f} {cy:.2f} {p1[0]:.2f} {p1[1]:.2f}" fill="none" stroke="{html.escape(rc, quote=True)}" stroke-width="{sw:.2f}" stroke-opacity="{r_alpha:.2f}"{dash_attr}/>'
                            )
                        else:
                            mx = (p0[0] + p1[0]) / 2.0 + 30.0
                            my = (p0[1] + p1[1]) / 2.0
                            parts.append(
                                f'<path d="M {p0[0]:.2f} {p0[1]:.2f} Q {mx:.2f} {my:.2f} {p1[0]:.2f} {p1[1]:.2f}" fill="none" stroke="{html.escape(rc, quote=True)}" stroke-width="{sw:.2f}" stroke-opacity="{r_alpha:.2f}"{dash_attr}/>'
                            )
                    elif gt == "tanglelink":
                        midx = (p0[0] + p1[0]) / 2.0
                        parts.append(
                            f'<path d="M {p0[0]:.2f} {p0[1]:.2f} C {midx:.2f} {p0[1]:.2f} {midx:.2f} {p1[1]:.2f} {p1[0]:.2f} {p1[1]:.2f}" fill="none" stroke="{html.escape(rc, quote=True)}" stroke-width="{sw:.2f}" stroke-opacity="{min(r_alpha, 0.75):.2f}"{dash_attr}/>'
                        )
                    elif gt == "range":
                        parts.append(
                            f'<line x1="{p0[0]:.2f}" y1="{p0[1]:.2f}" x2="{p1[0]:.2f}" y2="{p1[1]:.2f}" stroke="{html.escape(rc, quote=True)}" stroke-width="{max(sw, 3.0):.2f}" stroke-opacity="{r_alpha:.2f}" stroke-linecap="round"/>'
                        )
                    elif gt == "treescale":
                        parts.append(
                            f'<line x1="{p0[0]:.2f}" y1="{p0[1]:.2f}" x2="{p1[0]:.2f}" y2="{p1[1]:.2f}" stroke="{html.escape(rc, quote=True)}" stroke-width="{max(sw, 1.5):.2f}"/>'
                        )
                        parts.append(f'<line x1="{p0[0]:.2f}" y1="{p0[1] - 3:.2f}" x2="{p0[0]:.2f}" y2="{p0[1] + 3:.2f}" stroke="{html.escape(rc, quote=True)}" stroke-width="1.2"/>')
                        parts.append(f'<line x1="{p1[0]:.2f}" y1="{p1[1] - 3:.2f}" x2="{p1[0]:.2f}" y2="{p1[1] + 3:.2f}" stroke="{html.escape(rc, quote=True)}" stroke-width="1.2"/>')
                        lbl = r.get("label")
                        if not _is_scalar_na(lbl) and lbl is not None and not _is_scalar_na(r.get("label_x")) and not _is_scalar_na(r.get("label_y")):
                            lp = map_pt(float(r["label_x"]), float(r["label_y"]), panel=pan)
                            parts.append(
                                f'<text x="{lp[0]:.2f}" y="{lp[1] - 4:.2f}" text-anchor="middle" font-family="sans-serif" font-size="10" fill="{html.escape(rc, quote=True)}">{html.escape(str(lbl))}</text>'
                            )
                    else:
                        parts.append(
                            f'<line x1="{p0[0]:.2f}" y1="{p0[1]:.2f}" x2="{p1[0]:.2f}" y2="{p1[1]:.2f}" stroke="{html.escape(rc, quote=True)}" stroke-width="{sw:.2f}" stroke-opacity="{r_alpha:.2f}"{dash_attr}/>'
                        )

            # 5. Inset Charts (nodepie / nodebar)
            elif gt == "inset" and isinstance(layer.params.get("insets"), dict):
                insets_map: Dict[Any, Any] = layer.params["insets"]
                in_w = float(layer.params.get("width", 0.1) or 0.1)
                in_h = float(layer.params.get("height", 0.1) or 0.1)
                pie_r = max(9.0, min(26.0, min(plot_w, plot_h) * max(in_w, in_h) * 0.45))
                node_rows = {int(r["node"]): r for r in self.data.to_records() if not _is_scalar_na(r.get("node"))}
                x_col = "branch" if layer.params.get("x") == "branch" and "branch" in self.data else "x"
                for nid, spec in insets_map.items():
                    nr = node_rows.get(int(nid)) if isinstance(nid, (int, np.integer, str)) else None
                    if nr is None or _is_scalar_na(nr.get(x_col)) or _is_scalar_na(nr.get("y")):
                        continue
                    px, py = map_pt(float(nr[x_col]), float(nr["y"]))
                    vals_dict = spec.get("values", {})
                    cols_dict = spec.get("color") or {}
                    ialpha = float(spec.get("alpha", 0.85) or 0.85)
                    if spec.get("type") == "pie" and isinstance(vals_dict, dict):
                        total = sum(max(float(v), 0.0) for v in vals_dict.values())
                        if total <= 0:
                            continue
                        cur_a = -math.pi / 2.0
                        keys_list = list(vals_dict.keys())
                        for ki, k in enumerate(keys_list):
                            val = max(float(vals_dict[k]), 0.0)
                            frac = val / total
                            if frac <= 0:
                                continue
                            scol = _resolve_color_value(cols_dict.get(k), _DISCRETE_PALETTE[ki % len(_DISCRETE_PALETTE)])
                            if frac >= 0.9999:
                                parts.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="{pie_r:.2f}" fill="{html.escape(scol, quote=True)}" fill-opacity="{ialpha:.2f}" stroke="white" stroke-width="0.6"/>')
                            else:
                                nxt_a = cur_a + frac * 2.0 * math.pi
                                x0_p, y0_p = px + pie_r * math.cos(cur_a), py + pie_r * math.sin(cur_a)
                                x1_p, y1_p = px + pie_r * math.cos(nxt_a), py + pie_r * math.sin(nxt_a)
                                large_arc = 1 if frac > 0.5 else 0
                                parts.append(
                                    f'<path d="M {px:.2f} {py:.2f} L {x0_p:.2f} {y0_p:.2f} A {pie_r:.2f} {pie_r:.2f} 0 {large_arc} 1 {x1_p:.2f} {y1_p:.2f} Z" fill="{html.escape(scol, quote=True)}" fill-opacity="{ialpha:.2f}" stroke="white" stroke-width="0.6"/>'
                                )
                                cur_a = nxt_a
                    elif spec.get("type") == "bar" and isinstance(vals_dict, dict):
                        keys_list = list(vals_dict.keys())
                        pos_mode = str(spec.get("position", "stack") or "stack")
                        if pos_mode == "stack":
                            bw_tot = pie_r * 2.4
                            bh_tot = pie_r * 0.95
                            tot_v = sum(max(float(vals_dict[k]), 0.0) for k in keys_list) or 1.0
                            cur_x = px - bw_tot / 2.0
                            top_y = py - bh_tot / 2.0
                            for ki, k in enumerate(keys_list):
                                val = max(float(vals_dict[k]), 0.0)
                                seg_w = (val / tot_v) * bw_tot
                                if seg_w <= 0.1:
                                    continue
                                scol = _resolve_color_value(cols_dict.get(k), _DISCRETE_PALETTE[ki % len(_DISCRETE_PALETTE)])
                                parts.append(
                                    f'<rect x="{cur_x:.2f}" y="{top_y:.2f}" width="{seg_w:.2f}" height="{bh_tot:.2f}" fill="{html.escape(scol, quote=True)}" fill-opacity="{ialpha:.2f}" stroke="white" stroke-width="0.5"/>'
                                )
                                cur_x += seg_w
                        else:
                            bw_tot = pie_r * 2.0
                            bh_tot = pie_r * 1.6
                            mx_v = max((abs(float(vals_dict[k])) for k in keys_list), default=1.0) or 1.0
                            sub_w = bw_tot / max(len(keys_list), 1)
                            for ki, k in enumerate(keys_list):
                                val = max(float(vals_dict[k]), 0.0)
                                bh_i = (val / mx_v) * bh_tot
                                scol = _resolve_color_value(cols_dict.get(k), _DISCRETE_PALETTE[ki % len(_DISCRETE_PALETTE)])
                                rx = px - bw_tot / 2.0 + ki * sub_w
                                ry = py + bh_tot / 2.0 - bh_i
                                parts.append(
                                    f'<rect x="{rx:.2f}" y="{ry:.2f}" width="{max(sub_w - 1.0, 2.0):.2f}" height="{bh_i:.2f}" fill="{html.escape(scol, quote=True)}" fill-opacity="{ialpha:.2f}"/>'
                                )

            # 6. Points (point, point2, tippoint, nodepoint, rootpoint)
            elif gt in ("point", "point2", "tippoint", "nodepoint", "rootpoint") and len(ldf) > 0:
                recs = ldf.to_records()
                has_c_map = "color" in merged_aes or "colour" in merged_aes
                cmap_fn = None
                if has_c_map:
                    c_expr6 = merged_aes.get("color", merged_aes.get("colour"))
                    cvals = [r.get("color", r.get("colour")) for r in recs]
                    cmap_fn, leg_info = _build_value_color_mapper(
                        cvals,
                        layer.params.get("_color_scale") or self.scales.get("color"),
                        default_color=color,
                        force_discrete=(str(c_expr6) in ("group", "focus")),
                        levels=(ldf._levels.get("color") or ldf._levels.get("colour") or self.data._levels.get(str(c_expr6))),
                    )
                    _add_legend_unique(legends_to_draw, leg_info)
                r_default = max(size * (1.5 if compact else 1.8), 2.5 if compact else 3.2)
                has_s_map = "size" in merged_aes
                s_min = s_max = 0.0
                r_lo, r_hi = (2.2, 6.0)
                if has_s_map:
                    sz_scale = self.scales.get("size")
                    if sz_scale is not None and isinstance(sz_scale.params.get("range"), (list, tuple)) and len(sz_scale.params["range"]) >= 2:
                        r_lo, r_hi = float(sz_scale.params["range"][0]), float(sz_scale.params["range"][1])
                    svals = []
                    for r in recs:
                        sv = r.get("size")
                        if not _is_scalar_na(sv):
                            try:
                                svals.append(float(sv))
                            except Exception:
                                pass
                    if svals:
                        s_min, s_max = min(svals), max(svals)
                for r in recs:
                    pan = str(r.get("PANEL", default_pan) or default_pan)
                    if pan not in panel_boxes:
                        pan = default_pan
                    if _is_scalar_na(r.get("x")) or _is_scalar_na(r.get("y")):
                        continue
                    px, py = map_pt(float(r["x"]), float(r["y"]), panel=pan)
                    rc = cmap_fn(r.get("color", r.get("colour"))) if cmap_fn else _resolve_color_value(r.get("color", r.get("colour")), color)
                    if has_s_map and not _is_scalar_na(r.get("size")):
                        try:
                            sv = float(r["size"])
                            if s_max > s_min:
                                mapped_sz = r_lo + (sv - s_min) / (s_max - s_min) * (r_hi - r_lo)
                            else:
                                mapped_sz = (r_lo + r_hi) / 2.0
                            r_px = mapped_sz * (1.4 if compact else 1.8)
                        except Exception:
                            r_px = r_default
                    else:
                        r_sz = float(r.get("size", size) or size)
                        r_px = max(r_sz * (1.4 if compact else 1.8), 2.5) if "size" in r else r_default
                    parts.append(
                        f'<circle cx="{px:.2f}" cy="{py:.2f}" r="{r_px:.2f}" fill="{html.escape(rc, quote=True)}" fill-opacity="{alpha:.2f}" stroke="white" stroke-width="0.5"/>'
                    )

            # 7. Text & Labels
            elif gt in ("tiplab", "tiplab2", "nodelab", "nodelab2", "text", "text2", "label", "label2", "cladelab_text", "strip_text") and len(ldf) > 0:
                base_fsize = float(layer.params.get("fontsize", size if size != 1.1 else 3.2) or 3.2)
                min_f = 7 if (compact or len(tip_ys) >= 45) else 9
                fsize = max(int(round(base_fsize * (2.5 if compact else 3.1))), min_f)
                if not is_polar and not is_flipped and len(tip_ys) >= 40 and gt in ("tiplab", "tiplab2", "nodelab", "nodelab2"):
                    row_px = plot_h / max(len(tip_ys), 1)
                    fsize = max(6, min(fsize, int(row_px * 0.72)))
                is_box_label = gt in ("label", "label2") or layer.params.get("geom") == "label"
                recs7 = ldf.to_records()
                cmap7 = None
                if ("color" in merged_aes or "colour" in merged_aes) and recs7:
                    c_expr7 = merged_aes.get("color", merged_aes.get("colour"))
                    cvals7 = [r.get("color", r.get("colour")) for r in recs7]
                    cmap7, leg7 = _build_value_color_mapper(
                        cvals7,
                        layer.params.get("_color_scale") or self.scales.get("color"),
                        default_color=color,
                        force_discrete=(str(c_expr7) in ("group", "focus")),
                        levels=(ldf._levels.get("color") or ldf._levels.get("colour") or self.data._levels.get(str(c_expr7))),
                    )
                    _add_legend_unique(legends_to_draw, leg7)
                seen_tiplabs = set()
                for r in recs7:
                    pan = str(r.get("PANEL", default_pan) or default_pan)
                    if pan not in panel_boxes:
                        pan = default_pan
                    lbl_raw = r.get("label")
                    if _is_scalar_na(lbl_raw) or lbl_raw is None or str(lbl_raw) == "":
                        continue
                    if _is_scalar_na(r.get("x")) or _is_scalar_na(r.get("y")):
                        continue
                    parse_f = bool(r.get("parse", layer.params.get("parse", False)))
                    lbl_clean, f_style, f_weight = _parse_plotmath_label(str(lbl_raw), parse_f)
                    if not lbl_clean:
                        continue
                    rx_val, ry_val = float(r["x"]), float(r["y"])
                    px, py = map_pt(rx_val, ry_val, panel=pan)
                    if gt in ("tiplab", "tiplab2"):
                        tip_key = (pan, round(px, 1), round(py, 1), lbl_clean)
                        if tip_key in seen_tiplabs:
                            continue
                        seen_tiplabs.add(tip_key)
                    rc = cmap7(r.get("color", r.get("colour"))) if cmap7 else _resolve_color_value(r.get("color", r.get("colour")), color)
                    hj = float(r.get("hjust", layer.params.get("hjust", 0.0)) or 0.0)
                    style_attrs = ""
                    if f_style != "normal":
                        style_attrs += f' font-style="{f_style}"'
                    if f_weight != "normal":
                        style_attrs += f' font-weight="{f_weight}"'

                    if is_polar:
                        deg = polar_deg(ry_val) % 360.0
                        if gt in ("text", "text2") and ry_val < tip_y_lo - 0.15 and "angle" in layer.params:
                            col_fs = max(7, min(fsize, 9))
                            rot_col = -90.0 if (deg < 90.0 or deg > 270.0) else 90.0
                            esc_lbl = html.escape(lbl_clean)
                            parts.append(
                                f'<text x="{px:.2f}" y="{py + col_fs * 0.34:.2f}" text-anchor="middle" transform="rotate({rot_col:.1f},{px:.2f},{py:.2f})" font-family="sans-serif" font-size="{col_fs}" fill="white" stroke="white" stroke-width="2.5" stroke-linejoin="round" opacity="0.88"{style_attrs}>{esc_lbl}</text>'
                            )
                            parts.append(
                                f'<text x="{px:.2f}" y="{py + col_fs * 0.34:.2f}" text-anchor="middle" transform="rotate({rot_col:.1f},{px:.2f},{py:.2f})" font-family="sans-serif" font-size="{col_fs}" fill="{html.escape(rc, quote=True)}"{style_attrs}>{esc_lbl}</text>'
                            )
                            continue
                        if 90.0 < deg < 270.0:
                            rot_svg = 180.0 - deg
                            anchor = "end"
                            pad_r = -4.0
                        else:
                            rot_svg = -deg
                            anchor = "start"
                            pad_r = 4.0
                        th = math.radians(deg)
                        tx = px + 4.0 * math.cos(th)
                        ty = py - 4.0 * math.sin(th)
                        parts.append(
                            f'<text x="{tx:.2f}" y="{ty + fsize * 0.32:.2f}" text-anchor="{anchor}" transform="rotate({rot_svg:.2f},{tx:.2f},{ty:.2f})" font-family="sans-serif" font-size="{fsize}" fill="{html.escape(rc, quote=True)}"{style_attrs}>{html.escape(lbl_clean)}</text>'
                        )
                    elif is_flipped:
                        rot_svg = -90.0
                        parts.append(
                            f'<text x="{px:.2f}" y="{py + 5:.2f}" text-anchor="end" transform="rotate({rot_svg:.1f},{px:.2f},{py + 5:.2f})" font-family="sans-serif" font-size="{fsize}" fill="{html.escape(rc, quote=True)}"{style_attrs}>{html.escape(lbl_clean)}</text>'
                        )
                    else:
                        if "angle" in layer.mapping or gt in ("cladelab_text", "strip_text"):
                            ang = float(r.get("angle", layer.params.get("angle", 0.0)) or 0.0)
                        elif "angle" in layer.params and layer.params.get("angle") is not None and layer.params.get("angle") != "auto":
                            ang = float(layer.params["angle"])
                        else:
                            ang = 0.0
                        if hj >= 0.75:
                            anchor = "end"
                            tx = px - 4.0
                        elif hj >= 0.25:
                            anchor = "middle"
                            tx = px
                        else:
                            anchor = "start"
                            tx = px + 4.0
                        ty = py + fsize * 0.34
                        rot_attr = f' transform="rotate({-ang:.2f},{tx:.2f},{ty:.2f})"' if abs(ang) > 0.1 else ""
                        if is_box_label:
                            bw_pill = len(lbl_clean) * fsize * 0.62 + 8.0
                            bh_pill = fsize + 6.0
                            bx_pill = tx - (bw_pill / 2.0 if anchor == "middle" else (bw_pill if anchor == "end" else 4.0))
                            parts.append(
                                f'<rect x="{bx_pill:.2f}" y="{py - bh_pill / 2.0:.2f}" width="{bw_pill:.2f}" height="{bh_pill:.2f}" rx="2.5" fill="white" stroke="{html.escape(rc, quote=True)}" stroke-width="0.8"/>'
                            )
                        parts.append(
                            f'<text x="{tx:.2f}" y="{ty:.2f}" text-anchor="{anchor}"{rot_attr} font-family="sans-serif" font-size="{fsize}" fill="{html.escape(rc, quote=True)}"{style_attrs}>{html.escape(lbl_clean)}</text>'
                        )

        # Single-panel theme axes (theme_tree2 / theme_dendrogram)
        if not is_multi_panel and not is_polar and (show_x_axis or show_y_axis):
            bx, by, bw, bh, px0, px1 = panel_boxes["Tree"]
            if is_flipped:
                # Vertical height axis on the left for dendrogram
                parts.append(f'<line x1="{bx - 4:.2f}" y1="{by:.2f}" x2="{bx - 4:.2f}" y2="{by + bh:.2f}" stroke="#1F2937" stroke-width="1.2"/>')
                # Tree x goes from -span (top) to 0 (bottom); label height above 0 as positive
                for tv, _ in _nice_ticks(px0, px1, n=5):
                    _, ty = map_pt(tv, y_min)
                    h_val = abs(tv - px1)
                    hlbl = f"{h_val:g}"
                    parts.append(f'<line x1="{bx - 8:.2f}" y1="{ty:.2f}" x2="{bx - 4:.2f}" y2="{ty:.2f}" stroke="#1F2937" stroke-width="1"/>')
                    parts.append(f'<text x="{bx - 11:.2f}" y="{ty + 3:.2f}" text-anchor="end" font-family="sans-serif" font-size="9.5" fill="#374151">{html.escape(hlbl)}</text>')
            elif show_x_axis:
                ay = by + bh + 6.0
                parts.append(f'<line x1="{bx:.2f}" y1="{ay:.2f}" x2="{bx + bw:.2f}" y2="{ay:.2f}" stroke="#1F2937" stroke-width="1.2"/>')
                x_lbl_fn = self.scales["x"].params.get("labels") if "x" in self.scales else None
                for tv, tlbl in _nice_ticks(px0, px1, n=6, label_fn=x_lbl_fn):
                    tx, _ = map_pt(tv, y_min)
                    parts.append(f'<line x1="{tx:.2f}" y1="{ay:.2f}" x2="{tx:.2f}" y2="{ay + 4:.2f}" stroke="#1F2937" stroke-width="1"/>')
                    parts.append(f'<text x="{tx:.2f}" y="{ay + 15:.2f}" text-anchor="middle" font-family="sans-serif" font-size="9.5" fill="#374151">{html.escape(tlbl)}</text>')

        # Compact right-side stacked legends for heatmaps / MSA when not in tight grid cell
        if legends_to_draw and not compact:
            lx = w - 100
            ly = max(pad_top + 6, 44)
            for leg_idx, leg in enumerate(legends_to_draw):
                if ly > h - 40:
                    break
                ltitle = str(leg.get("title", "") or "value")
                parts.append(f'<text x="{lx:.1f}" y="{ly:.1f}" font-family="sans-serif" font-size="10.5" font-weight="bold" fill="#1F2937">{html.escape(ltitle)}</text>')
                ly += 10
                if leg.get("type") == "discrete":
                    for item_lbl, item_col in list(leg.get("items", []))[:10]:
                        parts.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="11" height="11" fill="{html.escape(item_col, quote=True)}" stroke="#9CA3AF" stroke-width="0.4"/>')
                        parts.append(f'<text x="{lx + 15:.1f}" y="{ly + 9.5:.1f}" font-family="sans-serif" font-size="9.5" fill="#374151">{html.escape(str(item_lbl)[:14])}</text>')
                        ly += 14
                    ly += 18
                elif leg.get("type") == "continuous":
                    stops = leg.get("stops", ["#440154", "#fde725"])
                    gid = f"gggrad_{leg_idx}_{abs(hash(tuple(stops))) % 100000}"
                    stop_tags = "".join(
                        f'<stop offset="{(i / max(len(stops) - 1, 1)) * 100:.1f}%" stop-color="{html.escape(sc, quote=True)}"/>'
                        for i, sc in enumerate(stops)
                    )
                    parts.append(f'<defs><linearGradient id="{gid}" x1="0" y1="1" x2="0" y2="0">{stop_tags}</linearGradient></defs>')
                    parts.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="12" height="72" fill="url(#{gid})" stroke="#9CA3AF" stroke-width="0.5"/>')
                    parts.append(f'<text x="{lx + 16:.1f}" y="{ly + 8:.1f}" font-family="sans-serif" font-size="9" fill="#374151">{leg.get("vmax", 1.0):.2g}</text>')
                    parts.append(f'<text x="{lx + 16:.1f}" y="{ly + 71:.1f}" font-family="sans-serif" font-size="9" fill="#374151">{leg.get("vmin", 0.0):.2g}</text>')
                    ly += 88

        if _interactive_overlay:
            append_interactive_node_overlay(
                parts,
                self.data,
                map_pt,
                layout_str=str(self.layout),
                is_polar=is_polar,
                polar_arc_d_fn=polar_arc_d if is_polar else None,
            )
        parts.append("</svg>")
        return "\n".join(parts)

    def _repr_svg_(self) -> str:
        return self.to_svg()

    def _render_rgba(self, width: int = 800, height: int = 600) -> np.ndarray:
        """Fallback pure-numpy RGBA rasterizer when rsvg-convert is unavailable."""
        w = max(int(width), 1)
        h = max(int(height), 1)
        img = np.full((h, w, 4), 255, dtype=np.uint8)

        built = self.build()
        pad_left = min(40, w // 6)
        pad_right = min(60, w // 5)
        pad_top = min(30, h // 6)
        pad_bottom = min(35, h // 6)
        plot_w = max(w - pad_left - pad_right, 10)
        plot_h = max(h - pad_top - pad_bottom, 10)
        x_min, x_max = built.layout.panel_scales_x[0].range.range
        y_min, y_max = built.layout.panel_scales_y[0].range.range
        dx = (x_max - x_min) if x_max != x_min else 1.0
        dy = (y_max - y_min) if y_max != y_min else 1.0

        def map_pt(x: float, y: float) -> Tuple[int, int]:
            px = int(round(pad_left + ((x - x_min) / dx) * plot_w))
            py = int(round(pad_top + plot_h - ((y - y_min) / dy) * plot_h))
            return (max(0, min(w - 1, px)), max(0, min(h - 1, py)))

        def draw_line(x0: int, y0: int, x1: int, y1: int, rgba: Tuple[int, int, int, int]) -> None:
            steps = max(abs(x1 - x0), abs(y1 - y0), 1)
            xs = np.linspace(x0, x1, steps + 1).astype(int)
            ys = np.linspace(y0, y1, steps + 1).astype(int)
            valid = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
            img[ys[valid], xs[valid]] = rgba

        recs = self.data.to_records()
        node_map = {r.get("node"): r for r in recs if r.get("node") is not None}
        for r in recs:
            nid, pid = r.get("node"), r.get("parent")
            if nid is None or pid is None or nid == pid or pid not in node_map:
                continue
            prow = node_map[pid]
            if _is_scalar_na(r.get("x")) or _is_scalar_na(r.get("y")) or _is_scalar_na(prow.get("x")) or _is_scalar_na(prow.get("y")):
                continue
            p_parent = map_pt(float(prow["x"]), float(prow["y"]))
            p_corner = map_pt(float(prow["x"]), float(r["y"]))
            p_child = map_pt(float(r["x"]), float(r["y"]))
            draw_line(p_parent[0], p_parent[1], p_corner[0], p_corner[1], (30, 30, 30, 255))
            draw_line(p_corner[0], p_corner[1], p_child[0], p_child[1], (30, 30, 30, 255))

        return img

    def to_png(self, width: int = 800, height: int = 600) -> bytes:
        """Render the GGTree plot to valid 8-bit RGBA PNG bytes."""
        svg = self.to_svg(width=width, height=height)
        png = _svg_to_png_bytes(svg, int(width), int(height))
        if png is not None:
            return png
        return _encode_png_rgba(self._render_rgba(width=width, height=height))

    def to_pdf(
        self,
        width: float = 800,
        height: float = 600,
        page_width_pt: Optional[float] = None,
        page_height_pt: Optional[float] = None,
    ) -> bytes:
        """Render the plot as a genuine multi-layer vector PDF document."""
        w = max(int(round(width)), 100)
        h = max(int(round(height)), 80)
        svg = self.to_svg(width=w, height=h)
        return svg_to_pdf_bytes(svg, width=w, height=h, page_width_pt=page_width_pt, page_height_pt=page_height_pt)

    def to_html(
        self,
        width: int = 860,
        height: int = 620,
        interactive: bool = True,
        full_html: bool = True,
        title: Optional[str] = None,
    ) -> str:
        """Render the plot as self-contained interactive (or static) HTML."""
        return render_interactive_html(
            self,
            width=width,
            height=height,
            interactive=interactive,
            full_html=full_html,
            title=title,
        )

    def to_interactive(
        self,
        width: int = 860,
        height: int = 620,
        title: Optional[str] = None,
    ) -> InteractivePlot:
        """Wrap the plot in a Jupyter-displayable InteractivePlot container."""
        return InteractivePlot(self, width=width, height=height, title=title)

    def save(self, filename: Union[str, Path], width: float = 8, height: float = 6, dpi: int = 100, **kwargs: Any) -> str:
        return ggsave(filename, plot=self, width=width, height=height, dpi=dpi, **kwargs)

    def __repr__(self) -> str:
        return f"<GGTree layout={self.layout!r} nodes={len(self.data)} layers={len(self.layers)}>"


# ---------------------------------------------------------------------------
# Color Scale, Label Parsing, Axis Ticks, and Vector/Raster Rendering Helpers
# ---------------------------------------------------------------------------
_VIRIDIS_PALETTES: Dict[str, List[str]] = {
    "D": ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"],
    "A": ["#000004", "#3b0f70", "#8c2981", "#de4968", "#fe9f6d", "#fcfdbf"],
    "B": ["#000004", "#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4"],
    "C": ["#0d0887", "#6a00a8", "#b12a90", "#e16462", "#fca636", "#f0f921"],
    "E": ["#00224e", "#4b5263", "#7c7b78", "#b5a86d", "#fee838"],
}

_DISCRETE_PALETTE: List[str] = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    "#7E6148", "#B09C85", "#1F78B4", "#33A02C",
    "#FF7F00", "#6A3D9A", "#E31A1C", "#B15928",
]

_BREWER_PALETTES: Dict[str, List[str]] = {
    "Set1": ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999"],
    "Set2": ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"],
    "Dark2": ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"],
    "Paired": ["#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00"],
}

_MSA_PALETTE: Dict[str, str] = {
    "A": "#4DAF4A", "C": "#377EB8", "G": "#FF7F00", "T": "#E41A1C", "U": "#E41A1C",
    "R": "#66C2A5", "N": "#8DA0CB", "D": "#E78AC3", "E": "#A6D854", "Q": "#FFD92F",
    "H": "#E5C494", "I": "#1B9E77", "L": "#D95F02", "K": "#7570B3", "M": "#E7298A",
    "F": "#66A61E", "P": "#E6AB02", "S": "#A6761D", "W": "#666666", "Y": "#386CB0",
    "V": "#F0027F",
}

_KNOWN_CSS_COLORS = {
    "red", "green", "blue", "black", "white", "grey", "gray", "grey50", "gray50",
    "grey70", "gray70", "steelblue", "firebrick", "darkgreen", "darkred", "darkblue",
    "orange", "purple", "yellow", "pink", "cyan", "magenta", "brown", "navy",
    "lightblue", "light blue", "gold", "salmon", "tomato", "coral", "teal",
    "indigo", "violet", "olive", "tan", "beige", "none", "transparent", "skyblue",
    "forestgreen", "limegreen", "darkorange", "dodgerblue", "royalblue", "slategray",
}


def _is_direct_css_color(val: Any) -> bool:
    if _is_scalar_na(val) or val is None:
        return False
    s = str(val).strip()
    if not s:
        return False
    return s.startswith("#") or s.startswith("rgb(") or s.startswith("rgba(") or s.lower() in _KNOWN_CSS_COLORS


def _resolve_color_value(val: Any, default_color: str = "#222222") -> str:
    if _is_scalar_na(val) or val is None:
        return default_color
    s = str(val).strip()
    if not s:
        return default_color
    if _is_direct_css_color(s):
        return s
    return default_color


def _add_legend_unique(legends_list: List[Dict[str, Any]], leg_info: Optional[Dict[str, Any]]) -> None:
    if not leg_info:
        return
    key = (
        leg_info.get("type"),
        leg_info.get("title"),
        tuple(leg_info.get("items", ())) if leg_info.get("type") == "discrete" else (leg_info.get("vmin"), leg_info.get("vmax")),
    )
    for existing in legends_list:
        ekey = (
            existing.get("type"),
            existing.get("title"),
            tuple(existing.get("items", ())) if existing.get("type") == "discrete" else (existing.get("vmin"), existing.get("vmax")),
        )
        if ekey == key:
            return
    legends_list.append(leg_info)


def _hex_to_rgb_tuple(hex_color: str) -> Tuple[int, int, int]:
    from .themes_scales import _hex_to_rgb
    return _hex_to_rgb(hex_color)


def _interp_stops(t: float, stops: Sequence[str]) -> str:
    if not stops:
        return "#4682B4"
    if len(stops) == 1:
        return stops[0]
    t = max(0.0, min(1.0, float(t)))
    pos = t * (len(stops) - 1)
    idx = min(int(math.floor(pos)), len(stops) - 2)
    frac = pos - idx
    r1, g1, b1 = _hex_to_rgb_tuple(stops[idx])
    r2, g2, b2 = _hex_to_rgb_tuple(stops[idx + 1])
    r = int(round(r1 + frac * (r2 - r1)))
    g = int(round(g1 + frac * (g2 - g1)))
    b = int(round(b1 + frac * (b2 - b1)))
    return f"#{r:02x}{g:02x}{b:02x}"


def _build_value_color_mapper(
    values: Sequence[Any],
    scale_obj: Any = None,
    default_color: str = "#4682B4",
    is_msa: bool = False,
    force_discrete: bool = False,
    levels: Optional[Sequence[Any]] = None,
) -> Tuple[Callable[[Any], str], Optional[Dict[str, Any]]]:
    non_na = [v for v in values if not _is_scalar_na(v) and str(v) != ""]
    if not non_na:
        return (lambda v: "#EEEEEE"), None

    st = getattr(scale_obj, "scale_type", None)
    sp = getattr(scale_obj, "params", {}) if scale_obj is not None else {}
    legend_title = sp.get("name") or ""

    if scale_obj is None and not is_msa and all(_is_direct_css_color(v) for v in non_na):
        return (lambda v: _resolve_color_value(v, default_color)), None

    is_numeric = all(isinstance(v, (int, float, np.number)) and not isinstance(v, bool) for v in non_na)
    if is_numeric and not force_discrete and st not in ("manual", "brewer", "viridis_d"):
        fvals = [float(v) for v in non_na if math.isfinite(float(v))]
        vmin = min(fvals) if fvals else 0.0
        vmax = max(fvals) if fvals else 1.0
        if st == "viridis_c":
            opt = str(sp.get("option", "D") or "D").upper()
            stops = _VIRIDIS_PALETTES.get(opt, _VIRIDIS_PALETTES["D"])
        elif st == "gradient2":
            stops = [str(sp.get("low", "blue")), str(sp.get("mid", "white")), str(sp.get("high", "red"))]
        elif st in ("gradient", "continuous"):
            stops = [str(sp.get("low", "#132B43")), str(sp.get("high", "#56B1F7"))]
        else:
            stops = _VIRIDIS_PALETTES["D"]
        span = (vmax - vmin) if vmax > vmin else 1.0

        def map_num(v: Any) -> str:
            if _is_scalar_na(v) or v is None:
                return "#EEEEEE"
            try:
                fv = float(v)
                if not math.isfinite(fv):
                    return "#EEEEEE"
            except Exception:
                return "#EEEEEE"
            return _interp_stops((fv - vmin) / span, stops)

        return map_num, {"type": "continuous", "title": legend_title or "value", "stops": stops, "vmin": vmin, "vmax": vmax}

    seen_vals: List[str] = []
    seen_set = set()
    for v in non_na:
        sv = str(v)
        if sv not in seen_set:
            seen_set.add(sv)
            seen_vals.append(sv)
    if levels:
        uniq: List[str] = []
        uniq_set = set()
        for lv in levels:
            slv = str(lv)
            if slv in seen_set and slv not in uniq_set:
                uniq_set.add(slv)
                uniq.append(slv)
        for sv in seen_vals:
            if sv not in uniq_set:
                uniq_set.add(sv)
                uniq.append(sv)
    else:
        uniq = sorted(seen_vals)

    cat_map: Dict[str, str] = {}
    if is_msa:
        for i, u in enumerate(uniq):
            cat_map[u] = _MSA_PALETTE.get(u.upper(), _DISCRETE_PALETTE[i % len(_DISCRETE_PALETTE)])
    elif st == "manual" and sp.get("values") is not None:
        man = sp["values"]
        if isinstance(man, dict):
            for i, u in enumerate(uniq):
                if u in man:
                    cat_map[u] = str(man[u])
                elif u.lstrip("-").isdigit() and int(u) in man:
                    cat_map[u] = str(man[int(u)])
                else:
                    cat_map[u] = _DISCRETE_PALETTE[i % len(_DISCRETE_PALETTE)]
        else:
            man_l = list(man)
            for i, u in enumerate(uniq):
                cat_map[u] = str(man_l[i % len(man_l)]) if man_l else _DISCRETE_PALETTE[i % len(_DISCRETE_PALETTE)]
    elif st == "brewer":
        pname = str(sp.get("palette", "Set1") or "Set1")
        pal = _BREWER_PALETTES.get(pname, _DISCRETE_PALETTE)
        for i, u in enumerate(uniq):
            cat_map[u] = pal[i % len(pal)]
    elif st == "viridis_d":
        opt = str(sp.get("option", "D") or "D").upper()
        stops = _VIRIDIS_PALETTES.get(opt, _VIRIDIS_PALETTES["D"])
        k = max(len(uniq) - 1, 1)
        for i, u in enumerate(uniq):
            cat_map[u] = _interp_stops(i / k, stops)
    else:
        for i, u in enumerate(uniq):
            cat_map[u] = _DISCRETE_PALETTE[i % len(_DISCRETE_PALETTE)]

    def map_cat(v: Any) -> str:
        if _is_scalar_na(v) or v is None or str(v) == "":
            return "#EEEEEE"
        return cat_map.get(str(v), default_color)

    leg_items = [(u, cat_map[u]) for u in uniq]
    return map_cat, {"type": "discrete", "title": legend_title or ("base" if is_msa else "group"), "items": leg_items}


def _parse_plotmath_label(text: str, parse_flag: bool) -> Tuple[str, str, str]:
    """Convert basic R plotmath wrappers like italic('Clade IL') into clean text + SVG font style."""
    if not parse_flag or not text:
        return text, "normal", "normal"
    s = text.strip()
    m_bi = re.match(r"^bolditalic\(\s*['\"](.*)['\"]\s*\)$", s)
    if m_bi:
        return m_bi.group(1), "italic", "bold"
    m_it = re.match(r"^italic\(\s*['\"](.*)['\"]\s*\)$", s)
    if m_it:
        return m_it.group(1), "italic", "normal"
    m_bf = re.match(r"^bold\(\s*['\"](.*)['\"]\s*\)$", s)
    if m_bf:
        return m_bf.group(1), "normal", "bold"
    return s, "normal", "normal"


def _nice_ticks(vmin: float, vmax: float, n: int = 5, label_fn: Any = None) -> List[Tuple[float, str]]:
    """Compute human-friendly axis tick positions and formatted labels inside [vmin, vmax]."""
    if not (math.isfinite(vmin) and math.isfinite(vmax)) or vmax <= vmin:
        return []
    span = vmax - vmin
    raw_step = span / max(n - 1, 1)
    mag = 10.0 ** math.floor(math.log10(raw_step))
    norm = raw_step / mag
    if norm <= 1.5:
        step = 1.0 * mag
    elif norm <= 3.0:
        step = 2.0 * mag
    elif norm <= 7.0:
        step = 5.0 * mag
    else:
        step = 10.0 * mag

    start = math.ceil((vmin - 1e-9) / step) * step
    ticks: List[Tuple[float, str]] = []
    cur = start
    while cur <= vmax + 1e-9 and len(ticks) < 12:
        val = 0.0 if abs(cur) < 1e-12 else cur
        if callable(label_fn):
            try:
                lbl = str(label_fn(val))
            except Exception:
                lbl = f"{val:g}"
        else:
            lbl = f"{val:.4g}" if abs(val) < 1000 and abs(val) >= 0.01 or val == 0 else f"{val:g}"
        ticks.append((val, lbl))
        cur += step
    return ticks


def _svg_to_png_bytes(svg: str, width: int, height: int) -> Optional[bytes]:
    """Rasterize SVG string to high-resolution PNG bytes via rsvg-convert when installed."""
    rsvg_bin = shutil.which("rsvg-convert")
    if rsvg_bin is None and Path("/opt/homebrew/bin/rsvg-convert").exists():
        rsvg_bin = "/opt/homebrew/bin/rsvg-convert"
    if rsvg_bin is None:
        return None
    try:
        proc = subprocess.run(
            [rsvg_bin, "-w", str(max(int(width), 1)), "-h", str(max(int(height), 1))],
            input=svg.encode("utf-8"),
            capture_output=True,
            check=False,
        )
        if proc.returncode == 0 and proc.stdout.startswith(b"\x89PNG\r\n\x1a\n"):
            return proc.stdout
    except Exception:
        return None
    return None


def _encode_png_rgba(img: np.ndarray) -> bytes:
    """Encode an (H, W, 4) uint8 RGBA image array as valid PNG bytes."""
    h, w = int(img.shape[0]), int(img.shape[1])
    raw_rows = b"".join(b"\x00" + img[y].tobytes() for y in range(h))

    def png_chunk(chunk_type: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(chunk_type + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + chunk_type + data + struct.pack(">I", crc)

    ihdr = struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0)
    idat = zlib.compress(raw_rows, level=6)
    return b"\x89PNG\r\n\x1a\n" + png_chunk(b"IHDR", ihdr) + png_chunk(b"IDAT", idat) + png_chunk(b"IEND", b"")


def render_plotlist_svg(
    plotlist: Any,
    width: int = 800,
    height: int = 600,
    _interactive_overlay: bool = False,
) -> str:
    """Render a PlotList grid of GGTree / nested PlotList plots into a combined SVG."""
    w = int(width)
    h = int(height)
    ncol = max(int(getattr(plotlist, "ncol", 1) or 1), 1)
    nrow = max(int(getattr(plotlist, "nrow", 1) or 1), 1)
    cell_w = max(w // ncol, 80)
    cell_h = max(h // nrow, 60)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        f'<rect width="{w}" height="{h}" fill="white"/>',
    ]
    labels = getattr(plotlist, "labels", None)
    for idx, p in enumerate(plotlist.plots):
        r = idx // ncol
        c = idx % ncol
        ox = c * cell_w
        oy = r * cell_h
        if hasattr(p, "to_svg"):
            if isinstance(p, GGTree):
                sub_svg = p.to_svg(
                    width=cell_w,
                    height=cell_h,
                    _in_grid=(cell_w < 500 or cell_h < 400),
                    _interactive_overlay=_interactive_overlay,
                )
            else:
                try:
                    sub_svg = p.to_svg(
                        width=cell_w,
                        height=cell_h,
                        _interactive_overlay=_interactive_overlay,
                    )
                except TypeError:
                    sub_svg = p.to_svg(width=cell_w, height=cell_h)
            inner = sub_svg.split(">", 1)[1].rsplit("</svg>", 1)[0] if "</svg>" in sub_svg else sub_svg
            parts.append(f'<g transform="translate({ox},{oy})">{inner}')
            if labels and idx < len(labels):
                parts.append(
                    f'<text x="10" y="18" font-family="sans-serif" font-size="13" font-weight="bold" fill="#111827">{html.escape(str(labels[idx]))}</text>'
                )
            parts.append("</g>")
    parts.append("</svg>")
    return "\n".join(parts)


def render_plotlist_rgba(plotlist: Any, width: int = 800, height: int = 600) -> np.ndarray:
    """Render a PlotList grid of plots into an (H, W, 4) uint8 RGBA array."""
    w = max(int(width), 1)
    h = max(int(height), 1)
    ncol = max(int(getattr(plotlist, "ncol", 1) or 1), 1)
    nrow = max(int(getattr(plotlist, "nrow", 1) or 1), 1)
    cell_w = max(w // ncol, 1)
    cell_h = max(h // nrow, 1)

    canvas = np.full((h, w, 4), 255, dtype=np.uint8)
    for idx, p in enumerate(plotlist.plots):
        r = idx // ncol
        c = idx % ncol
        if r >= nrow or c >= ncol:
            continue
        ox = c * cell_w
        oy = r * cell_h
        cw = (w - ox) if (c == ncol - 1) else min(cell_w, w - ox)
        ch = (h - oy) if (r == nrow - 1) else min(cell_h, h - oy)
        if cw <= 0 or ch <= 0:
            continue
        if hasattr(p, "_render_rgba"):
            sub = p._render_rgba(width=cw, height=ch)
            sh = min(ch, sub.shape[0])
            sw = min(cw, sub.shape[1])
            canvas[oy : oy + sh, ox : ox + sw] = sub[:sh, :sw]
    return canvas


def render_plotlist_png(plotlist: Any, width: int = 800, height: int = 600) -> bytes:
    """Render a PlotList grid of plots to valid 8-bit RGBA PNG bytes."""
    svg = render_plotlist_svg(plotlist, width=width, height=height)
    png = _svg_to_png_bytes(svg, int(width), int(height))
    if png is not None:
        return png
    return _encode_png_rgba(render_plotlist_rgba(plotlist, width=width, height=height))


# ---------------------------------------------------------------------------
# Top-level constructors & functions: ggplot, ggtree, ggplot_build, is_ggplot, ggsave
# ---------------------------------------------------------------------------

def ggplot(data: Any = None, mapping: Optional[Aes] = None, **kwargs: Any) -> GGTree:
    """Create a new ggplot / GGTree object."""
    if isinstance(data, Aes) and mapping is None:
        mapping = data
        data = None
    if data is None:
        df = DataFrame()
    elif isinstance(data, DataFrame):
        df = data.copy()
    elif isinstance(data, (Phylo, TreeData, MultiPhylo, TreeDataList)) or hasattr(data, "is_gtangle"):
        df = fortify(data, **kwargs)
    else:
        df = DataFrame(data)
    return GGTree(data=df, mapping=mapping, tr=data)


def ggtree(
    tr: Any,
    mapping: Optional[Aes] = None,
    layout: str = "rectangular",
    open_angle: float = 0.0,
    mrsd: Any = None,
    as_Date: bool = False,
    yscale: str = "none",
    yscale_mapping: Optional[Dict[Any, float]] = None,
    ladderize: bool = True,
    right: bool = False,
    branch_length: str = "branch.length",
    root_position: float = 0.0,
    xlim: Optional[Any] = None,
    **kwargs: Any,
) -> GGTree:
    """
    Main ggtree constructor mirroring R's `ggtree(...)`.
    Supports all 15+ layouts:
      'rectangular', 'roundrect', 'slanted', 'ellipse', 'dendrogram',
      'circular', 'fan', 'inward_circular', 'radial', 'equal_angle',
      'daylight', 'ape', 'unrooted', 'tree_and_leaf', 'tidy', or a custom callable.
    """
    branch_length = kwargs.pop("branch.length", branch_length)
    open_angle = float(kwargs.pop("open.angle", open_angle) or 0.0)
    as_Date = bool(kwargs.pop("as.Date", as_Date))
    yscale_mapping = kwargs.pop("yscale.mapping", yscale_mapping)
    root_position = float(kwargs.pop("root.position", root_position) or 0.0)

    if yscale != "none" and layout == "slanted":
        raise ValueError("slanted layout does not support yscale")

    extra_params = {
        "ladderize": ladderize,
        "right": right,
        "branch_length": branch_length,
        "mrsd": mrsd,
        "as_Date": as_Date,
        "yscale": yscale,
        "yscale_mapping": yscale_mapping,
        "root_position": root_position,
    }

    fortify_kw_names = {"optimize", "optimize_side", "gap", "mirror", "MAX_COUNT", "hang"}
    fortify_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in fortify_kw_names}

    df = fortify(
        tr,
        layout=layout,
        ladderize=ladderize,
        right=right,
        branch_length=branch_length,
        mrsd=mrsd,
        as_Date=as_Date,
        yscale=yscale,
        yscale_mapping=yscale_mapping,
        root_position=root_position,
        **fortify_kwargs,
    )

    if isinstance(df, DataFrame) and df.attrs.get("layout") in ("tanglegram", "external"):
        layout_str = str(df.attrs["layout"])
    elif callable(layout):
        layout_str = "external"
    else:
        layout_str = str(layout)
        if layout_str == "unrooted":
            layout_str = "daylight"

    is_multi = isinstance(tr, (MultiPhylo, TreeDataList)) or (isinstance(df, DataFrame) and ".id" in df.columns)
    base_mapping = Aes(x="x", y="y")
    if mapping is not None:
        base_mapping = modifyList(base_mapping, mapping)

    p = GGTree(data=df, mapping=base_mapping, layout=layout_str, tr=tr, extra_params=extra_params)
    p.attrs["layout"] = layout_str
    p.data.attrs["layout"] = layout_str
    p = p + geom_tree(layout=layout_str, multiPhylo=is_multi, **kwargs)
    p = p + theme_tree()

    if layout_str in ("circular", "radial"):
        p = p + coord_polar(theta="y", start=-math.pi / 2.0, direction=-1, clip="off")
    elif layout_str == "fan":
        from .manipulate import open_tree
        p = open_tree(p, angle=open_angle)
    elif layout_str == "inward_circular":
        from .manipulate import layout_inward_circular
        lim_arg = xlim[0] if isinstance(xlim, (list, tuple)) and len(xlim) == 1 else (xlim[1] if isinstance(xlim, (list, tuple)) and len(xlim) >= 2 else xlim)
        p = p + layout_inward_circular(xlim=lim_arg)
    elif layout_str == "dendrogram":
        from .manipulate import layout_dendrogram
        p = p + layout_dendrogram()

    set_last_plot(p)
    return p


def ggplot_build(plot: GGTree) -> BuiltPlot:
    """Build a GGTree plot for inspection and testing (compatible with ggplot2::ggplot_build)."""
    if not isinstance(plot, GGTree):
        raise TypeError(f"Expected GGTree, got {type(plot).__name__}")
    return plot.build()


def is_ggplot(x: Any) -> bool:
    """Return True if x is a GGTree / PlotList object."""
    from .compose import PlotList
    return isinstance(x, (GGTree, PlotList))


def ggsave(
    filename: Union[str, Path],
    plot: Optional[Any] = None,
    width: float = 8,
    height: float = 6,
    dpi: int = 100,
    device: Optional[str] = None,
    interactive: bool = True,
    **kwargs: Any,
) -> str:
    """Save a GGTree / PlotList / InteractivePlot to SVG, PNG, PDF, or interactive HTML."""
    if plot is None:
        plot = last_plot()
    if plot is None:
        raise ValueError("No plot provided and no last_plot() available.")
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    fmt = (str(device).strip().lstrip(".").lower() if device is not None else path.suffix.lower().lstrip("."))
    px_w = int(round(width * dpi))
    px_h = int(round(height * dpi))
    if fmt == "png":
        if not hasattr(plot, "to_png"):
            raise TypeError(f"Plot of type {type(plot).__name__} does not support PNG export.")
        path.write_bytes(plot.to_png(width=px_w, height=px_h))
    elif fmt == "pdf":
        svg_w = max(int(round(width * 100)), 100)
        svg_h = max(int(round(height * 100)), 80)
        pw_pt = float(width) * 72.0
        ph_pt = float(height) * 72.0
        if hasattr(plot, "to_pdf"):
            try:
                pdf_bytes = plot.to_pdf(width=svg_w, height=svg_h, page_width_pt=pw_pt, page_height_pt=ph_pt)
            except TypeError:
                pdf_bytes = plot.to_pdf(width=svg_w, height=svg_h)
        elif hasattr(plot, "to_svg"):
            pdf_bytes = svg_to_pdf_bytes(
                plot.to_svg(width=svg_w, height=svg_h),
                width=svg_w,
                height=svg_h,
                page_width_pt=pw_pt,
                page_height_pt=ph_pt,
            )
        else:
            raise TypeError(f"Plot of type {type(plot).__name__} does not support PDF export.")
        path.write_bytes(pdf_bytes)
    elif fmt in ("html", "htm"):
        html_str = render_interactive_html(
            plot,
            width=px_w,
            height=px_h,
            interactive=bool(interactive),
            full_html=True,
        )
        path.write_text(html_str, encoding="utf-8")
    elif fmt == "svg":
        svg = plot.to_svg(width=px_w, height=px_h)
        path.write_text(svg, encoding="utf-8")
    else:
        raise ValueError(
            f"Unsupported output format {fmt!r} for {str(filename)!r}. Supported formats: svg, png, pdf, html."
        )
    return str(path)
