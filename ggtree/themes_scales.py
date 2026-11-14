"""
Themes, scales, coordinates, facets, and axis utilities for ggtree:
theme_tree, theme_tree2, theme_dendrogram, theme_inset, theme,
xlim, ylim, xlim_tree, xlim_expand, hexpand, vexpand, ggexpand,
scale_x_range, revts, scale_x_ggtree, scale_color, add_colorbar,
scale_color_manual, scale_fill_manual, scale_fill_gradient, etc.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union
import numpy as np

from .dataframe import DataFrame, _is_scalar_na


class Rel(float):
    """ggplot2::rel(x) relative sizing wrapper."""

    def __new__(cls, x: float) -> "Rel":
        return float.__new__(cls, float(x))

    def __repr__(self) -> str:
        return f"rel({float(self):g})"


def rel(x: float) -> Rel:
    """Relative size multiplier (ggplot2::rel)."""
    return Rel(x)


class Unit:
    """grid::unit(x, units) wrapper."""

    def __init__(self, x: Union[float, Sequence[float]], units: str = "cm"):
        self.x = x
        self.units = units

    def __repr__(self) -> str:
        return f"unit({self.x!r}, {self.units!r})"


def unit(x: Union[float, Sequence[float]], units: str = "cm") -> Unit:
    return Unit(x, units=units)


def margin(t: float = 0, r: float = 0, b: float = 0, l: float = 0, unit: str = "pt") -> Dict[str, Any]:
    return {"t": t, "r": r, "b": b, "l": l, "unit": unit}


def arrow(
    angle: float = 30,
    length: Any = None,
    ends: str = "last",
    type: str = "open",
) -> Dict[str, Any]:
    return {"angle": angle, "length": length or unit(0.25, "inches"), "ends": ends, "type": type}


def element_blank() -> Dict[str, Any]:
    return {"element": "blank"}


def element_text(**kwargs: Any) -> Dict[str, Any]:
    return {"element": "text", **kwargs}


def element_line(**kwargs: Any) -> Dict[str, Any]:
    return {"element": "line", **kwargs}


def element_rect(**kwargs: Any) -> Dict[str, Any]:
    return {"element": "rect", **kwargs}


def guide_legend(**kwargs: Any) -> Dict[str, Any]:
    return {"guide": "legend", **kwargs}


def guide_colorbar(**kwargs: Any) -> Dict[str, Any]:
    return {"guide": "colorbar", **kwargs}


class PlotModifier:
    """Base class for objects that can be added to a GGTree plot via `+`."""

    def apply(self, plot: Any) -> Any:
        return plot


class Theme(PlotModifier):
    """ggplot2 / ggtree Theme specification."""

    def __init__(self, name: str = "custom", bgcolor: str = "white", fgcolor: str = "black", show_x_axis: bool = False, **kwargs: Any):
        self.name = name
        self.bgcolor = bgcolor
        self.fgcolor = fgcolor
        self.show_x_axis = show_x_axis
        self.params: Dict[str, Any] = dict(kwargs)

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        if self.name != "custom":
            out.theme_spec["name"] = self.name
            out.theme_spec["bgcolor"] = self.bgcolor
            out.theme_spec["fgcolor"] = self.fgcolor
            out.theme_spec["show_x_axis"] = self.show_x_axis
        out.theme_spec.update(self.params)
        return out


def theme(**kwargs: Any) -> Theme:
    return Theme(name="custom", **kwargs)


def theme_tree(bgcolor: str = "white", fgcolor: str = "black", **kwargs: Any) -> Theme:
    """Tree theme with blank axes and clean background (R/theme.R)."""
    return Theme(name="tree", bgcolor=bgcolor, fgcolor=fgcolor, show_x_axis=False, **kwargs)


def theme_tree2(bgcolor: str = "white", fgcolor: str = "black", **kwargs: Any) -> Theme:
    """Tree theme with bottom x-axis line and ticks visible (R/theme.R)."""
    return Theme(name="tree2", bgcolor=bgcolor, fgcolor=fgcolor, show_x_axis=True, **kwargs)


def theme_dendrogram(bgcolor: str = "white", fgcolor: str = "black", **kwargs: Any) -> Theme:
    """Dendrogram theme with y-axis / height axis visible (R/theme.R)."""
    return Theme(name="dendrogram", bgcolor=bgcolor, fgcolor=fgcolor, show_x_axis=True, show_y_axis=True, **kwargs)


def theme_inset(**kwargs: Any) -> Theme:
    """Transparent inset theme for nodebar / nodepie (R/theme.R)."""
    return Theme(name="inset", bgcolor="transparent", show_x_axis=False, **kwargs)


def theme_bw(**kwargs: Any) -> Theme:
    return Theme(name="bw", bgcolor="white", fgcolor="black", show_x_axis=True, **kwargs)


def theme_minimal(**kwargs: Any) -> Theme:
    return Theme(name="minimal", bgcolor="white", fgcolor="black", show_x_axis=True, **kwargs)


def theme_void(**kwargs: Any) -> Theme:
    return Theme(name="void", bgcolor="white", fgcolor="black", show_x_axis=False, **kwargs)


def theme_classic(**kwargs: Any) -> Theme:
    return Theme(name="classic", bgcolor="white", fgcolor="black", show_x_axis=True, **kwargs)


class Labs(PlotModifier):
    def __init__(self, **kwargs: Any):
        self.labels = dict(kwargs)

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        out.labels.update(self.labels)
        return out


def labs(**kwargs: Any) -> Labs:
    return Labs(**kwargs)


def ggtitle(label: str, subtitle: Optional[str] = None) -> Labs:
    return Labs(title=label, subtitle=subtitle)


def xlab(label: str) -> Labs:
    return Labs(x=label)


def ylab(label: str) -> Labs:
    return Labs(y=label)


def guides(**kwargs: Any) -> PlotModifier:
    class _Guides(PlotModifier):
        def apply(self, plot: Any) -> Any:
            out = plot.copy()
            out.guides.update(kwargs)
            return out
    return _Guides()


# Axis & Limits Modifiers
class AxisLimits(PlotModifier):
    def __init__(self, axis: str, limits: Tuple[Optional[float], Optional[float]], panel: Optional[str] = None):
        self.axis = axis
        self.limits = limits
        self.panel = panel

    def apply(self, plot: Any) -> Any:
        from .geoms import Layer
        from .aes import aes
        out = plot.copy()
        lo, hi = self.limits
        if self.panel is not None:
            out.panel_limits.setdefault(self.panel, {})[self.axis] = (lo, hi)
            # Also add dummy blank layer with PANEL column like xlim_expand in R/axis.R
            vals = [v for v in (lo, hi) if v is not None and not _is_scalar_na(v)]
            if vals:
                d = DataFrame({"x": vals, "y": [1.0] * len(vals), ".panel": [self.panel] * len(vals)})
                out.layers.append(Layer(geom_type="blank", data=d, mapping=aes(x="x", y="y"), params={"panel": self.panel}))
        else:
            if self.axis == "x":
                out.xlim = (lo, hi)
            else:
                out.ylim = (lo, hi)
        return out


def _parse_lim_args(*args: Any) -> Tuple[Optional[float], Optional[float]]:
    if len(args) == 1:
        a = args[0]
        if isinstance(a, (list, tuple, np.ndarray)):
            if len(a) == 1:
                return (0.0, float(a[0]))
            lo = None if _is_scalar_na(a[0]) else float(a[0])
            hi = None if _is_scalar_na(a[1]) else float(a[1])
            return (lo, hi)
        return (0.0, float(a))
    if len(args) >= 2:
        lo = None if _is_scalar_na(args[0]) else float(args[0])
        hi = None if _is_scalar_na(args[1]) else float(args[1])
        return (lo, hi)
    return (None, None)


def xlim(*args: Any) -> AxisLimits:
    """Set x-axis limits, e.g. xlim(0, 6), xlim(NA, 6), or xlim(c(NA, 6))."""
    return AxisLimits("x", _parse_lim_args(*args))


def ylim(*args: Any) -> AxisLimits:
    """Set y-axis limits, e.g. ylim(0, 10) or ylim(c(0, 10))."""
    return AxisLimits("y", _parse_lim_args(*args))


def xlim_expand(xlim_vals: Any, panel: str = "Tree") -> AxisLimits:
    """Set x-axis limits for a specific facet panel (exact port of R/axis.R xlim_expand)."""
    return AxisLimits("x", _parse_lim_args(xlim_vals), panel=str(panel))


def xlim_tree(xlim_vals: Any) -> AxisLimits:
    """Set x-axis limits specifically for the Tree panel (exact port of R/axis.R xlim_tree)."""
    return xlim_expand(xlim_vals, panel="Tree")


class ExpandModifier(PlotModifier):
    """Port of ggexpand / hexpand / vexpand in R/axis.R."""

    def __init__(self, ratio: float, direction: int = 1, side: str = "h"):
        self.ratio = float(ratio)
        self.direction = int(direction)
        self.side = str(side).lower()

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        df = out.data
        if self.side in ("h", "horizontal", "x"):
            xs = [float(v) for v in df["x"].tolist() if not _is_scalar_na(v)]
            if xs:
                mn, mx = min(xs), max(xs)
                span = max(mx - mn, 1e-6)
                if self.direction >= 0:
                    cur_hi = out.xlim[1] if out.xlim and out.xlim[1] is not None else mx
                    out.xlim = (out.xlim[0] if out.xlim else mn, cur_hi + span * self.ratio)
                else:
                    cur_lo = out.xlim[0] if out.xlim and out.xlim[0] is not None else mn
                    out.xlim = (cur_lo - span * self.ratio, out.xlim[1] if out.xlim else mx)
        else:
            ys = [float(v) for v in df["y"].tolist() if not _is_scalar_na(v)]
            if ys:
                mn, mx = min(ys), max(ys)
                span = max(mx - mn, 1e-6)
                if self.direction >= 0:
                    cur_hi = out.ylim[1] if out.ylim and out.ylim[1] is not None else mx
                    out.ylim = (out.ylim[0] if out.ylim else mn, cur_hi + span * self.ratio)
                else:
                    cur_lo = out.ylim[0] if out.ylim and out.ylim[0] is not None else mn
                    out.ylim = (cur_lo - span * self.ratio, out.ylim[1] if out.ylim else mx)
        return out


def hexpand(ratio: float, direction: int = 1) -> ExpandModifier:
    return ExpandModifier(ratio=ratio, direction=direction, side="h")


def vexpand(ratio: float, direction: int = 1) -> ExpandModifier:
    return ExpandModifier(ratio=ratio, direction=direction, side="v")


def ggexpand(ratio: float, direction: int = 1, side: str = "h") -> ExpandModifier:
    return ExpandModifier(ratio=ratio, direction=direction, side=side)


def revts(treeview: Any = None) -> Any:
    """
    Reverse time scale of a ggtree plot so tips sit at 0 and ancestors are negative
    (exact port of revts in R/axis.R). Supports both `revts(p)` and `p + revts()`.
    """
    if treeview is None:
        class _RevtsMod(PlotModifier):
            def apply(self, plot: Any) -> Any:
                return revts(plot)
        return _RevtsMod()
    out = treeview.copy()
    xs = [float(v) for v in out.data["x"].tolist() if not _is_scalar_na(v)]
    mx = max(xs) if xs else 0.0
    out.data["x"] = [float(v) - mx if not _is_scalar_na(v) else None for v in out.data["x"].tolist()]
    if "branch" in out.data:
        out.data["branch"] = [float(v) - mx if not _is_scalar_na(v) else None for v in out.data["branch"].tolist()]
    out.data.attrs["revts.done"] = True
    return out


class ScaleXRange(PlotModifier):
    """Port of scale_x_range in R/axis.R."""

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        out.theme_spec["scale_x_range"] = True
        return out


def scale_x_range() -> ScaleXRange:
    return ScaleXRange()


# Coordinate systems & Facets
class Coord(PlotModifier):
    def __init__(self, coord_type: str, **kwargs: Any):
        self.coord_type = coord_type
        self.params = dict(kwargs)

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        out.coord = self
        if "xlim" in self.params and self.params["xlim"] is not None:
            out.xlim = _parse_lim_args(self.params["xlim"])
        if "ylim" in self.params and self.params["ylim"] is not None:
            out.ylim = _parse_lim_args(self.params["ylim"])
        return out


def coord_cartesian(xlim: Any = None, ylim: Any = None, expand: bool = True, clip: str = "on", **kwargs: Any) -> Coord:
    return Coord("cartesian", xlim=xlim, ylim=ylim, expand=expand, clip=clip, **kwargs)


def coord_flip(xlim: Any = None, ylim: Any = None, expand: bool = True, **kwargs: Any) -> Coord:
    return Coord("flip", xlim=xlim, ylim=ylim, expand=expand, **kwargs)


def coord_polar(theta: str = "x", start: float = 0.0, direction: int = 1, clip: str = "on", **kwargs: Any) -> Coord:
    return Coord("polar", theta=theta, start=start, direction=direction, clip=clip, **kwargs)


def coord_fixed(ratio: float = 1.0, xlim: Any = None, ylim: Any = None, expand: bool = True, **kwargs: Any) -> Coord:
    return Coord("fixed", ratio=ratio, xlim=xlim, ylim=ylim, expand=expand, **kwargs)


class FacetSpec(PlotModifier):
    def __init__(self, facet_type: str, facets: Any = None, **kwargs: Any):
        self.facet_type = facet_type
        self.facets = facets
        self.params = dict(kwargs)

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        out.facet_spec = self
        return out


def facet_wrap(facets: Any = "~.id", nrow: Optional[int] = None, ncol: Optional[int] = None, scales: str = "fixed", **kwargs: Any) -> FacetSpec:
    return FacetSpec("wrap", facets=facets, nrow=nrow, ncol=ncol, scales=scales, **kwargs)


def facet_grid(rows: Any = None, cols: Any = None, scales: str = "fixed", space: str = "fixed", **kwargs: Any) -> FacetSpec:
    return FacetSpec("grid", rows=rows, cols=cols, scales=scales, space=space, **kwargs)


# Scales & Color mapping
class Scale(PlotModifier):
    def __init__(self, aesthetic: str, scale_type: str, **kwargs: Any):
        self.aesthetic = "color" if aesthetic == "colour" else aesthetic
        self.scale_type = scale_type
        self.params = dict(kwargs)

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        out.scales[self.aesthetic] = self
        if self.aesthetic in ("fill", "color"):
            freeze_key = f"_{self.aesthetic}_scale_frozen"
            scale_key = f"_{self.aesthetic}_scale"
            for lyr in out.layers:
                if not lyr.params.get(freeze_key, False):
                    lyr.params[scale_key] = self
        if self.aesthetic == "x" and "limits" in self.params and self.params["limits"] is not None:
            out.xlim = _parse_lim_args(self.params["limits"])
        if self.aesthetic == "y" and "limits" in self.params and self.params["limits"] is not None:
            out.ylim = _parse_lim_args(self.params["limits"])
        return out


class NewScale(PlotModifier):
    """Port of ggnewscale::new_scale_fill / new_scale_color."""

    def __init__(self, aesthetic: str):
        self.aesthetic = "color" if aesthetic == "colour" else aesthetic

    def apply(self, plot: Any) -> Any:
        out = plot.copy()
        # Archive current scale and freeze existing layers so subsequent scales only apply to new layers
        archived_key = f"{self.aesthetic}_prev_{len(out.scales)}"
        if self.aesthetic in out.scales:
            out.scales[archived_key] = out.scales.pop(self.aesthetic)
        freeze_key = f"_{self.aesthetic}_scale_frozen"
        for lyr in out.layers:
            lyr.params[freeze_key] = True
        return out


def new_scale_fill() -> NewScale:
    return NewScale("fill")


def new_scale_color() -> NewScale:
    return NewScale("color")


def new_scale_colour() -> NewScale:
    return NewScale("color")


def scale_color_manual(values: Any = None, breaks: Any = None, labels: Any = None, name: Optional[str] = None, **kwargs: Any) -> Scale:
    return Scale("color", "manual", values=values, breaks=breaks, labels=labels, name=name, **kwargs)


def scale_colour_manual(values: Any = None, breaks: Any = None, labels: Any = None, name: Optional[str] = None, **kwargs: Any) -> Scale:
    return scale_color_manual(values=values, breaks=breaks, labels=labels, name=name, **kwargs)


def scale_fill_manual(values: Any = None, breaks: Any = None, labels: Any = None, name: Optional[str] = None, **kwargs: Any) -> Scale:
    return Scale("fill", "manual", values=values, breaks=breaks, labels=labels, name=name, **kwargs)


def scale_color_continuous(low: str = "#132B43", high: str = "#56B1F7", **kwargs: Any) -> Scale:
    return Scale("color", "continuous", low=low, high=high, **kwargs)


def scale_colour_continuous(low: str = "#132B43", high: str = "#56B1F7", **kwargs: Any) -> Scale:
    return scale_color_continuous(low=low, high=high, **kwargs)


def scale_fill_continuous(low: str = "#132B43", high: str = "#56B1F7", **kwargs: Any) -> Scale:
    return Scale("fill", "continuous", low=low, high=high, **kwargs)


def scale_color_gradient(low: str = "#132B43", high: str = "#56B1F7", **kwargs: Any) -> Scale:
    return Scale("color", "gradient", low=low, high=high, **kwargs)


def scale_colour_gradient(low: str = "#132B43", high: str = "#56B1F7", **kwargs: Any) -> Scale:
    return scale_color_gradient(low=low, high=high, **kwargs)


def scale_fill_gradient(low: str = "#132B43", high: str = "#56B1F7", **kwargs: Any) -> Scale:
    return Scale("fill", "gradient", low=low, high=high, **kwargs)


def scale_color_gradient2(low: str = "blue", mid: str = "white", high: str = "red", midpoint: float = 0.0, **kwargs: Any) -> Scale:
    return Scale("color", "gradient2", low=low, mid=mid, high=high, midpoint=midpoint, **kwargs)


def scale_fill_gradient2(low: str = "blue", mid: str = "white", high: str = "red", midpoint: float = 0.0, **kwargs: Any) -> Scale:
    return Scale("fill", "gradient2", low=low, mid=mid, high=high, midpoint=midpoint, **kwargs)


def scale_color_viridis_c(**kwargs: Any) -> Scale:
    return Scale("color", "viridis_c", **kwargs)


def scale_colour_viridis_c(**kwargs: Any) -> Scale:
    return Scale("color", "viridis_c", **kwargs)


def scale_fill_viridis_c(**kwargs: Any) -> Scale:
    return Scale("fill", "viridis_c", **kwargs)


def scale_color_viridis_d(**kwargs: Any) -> Scale:
    return Scale("color", "viridis_d", **kwargs)


def scale_colour_viridis_d(**kwargs: Any) -> Scale:
    return Scale("color", "viridis_d", **kwargs)


def scale_fill_viridis_d(**kwargs: Any) -> Scale:
    return Scale("fill", "viridis_d", **kwargs)


def scale_color_brewer(palette: str = "Set1", **kwargs: Any) -> Scale:
    return Scale("color", "brewer", palette=palette, **kwargs)


def scale_fill_brewer(palette: str = "Set1", **kwargs: Any) -> Scale:
    return Scale("fill", "brewer", palette=palette, **kwargs)


def scale_x_continuous(breaks: Any = None, labels: Any = None, limits: Any = None, expand: Any = None, **kwargs: Any) -> Scale:
    return Scale("x", "continuous", breaks=breaks, labels=labels, limits=limits, expand=expand, **kwargs)


def scale_y_continuous(breaks: Any = None, labels: Any = None, limits: Any = None, expand: Any = None, **kwargs: Any) -> Scale:
    return Scale("y", "continuous", breaks=breaks, labels=labels, limits=limits, expand=expand, **kwargs)


def scale_x_date(**kwargs: Any) -> Scale:
    return Scale("x", "date", **kwargs)


def scale_x_ggtree(breaks: Any = None, labels: Any = None, **kwargs: Any) -> Scale:
    """Port of ggtree::scale_x_ggtree in R/scales.R."""
    return Scale("x", "ggtree", breaks=breaks, labels=labels, **kwargs)


def scale_size(range: Tuple[float, float] = (1.0, 6.0), **kwargs: Any) -> Scale:
    return Scale("size", "continuous", range=range, **kwargs)


def scale_size_continuous(range: Tuple[float, float] = (1.0, 6.0), **kwargs: Any) -> Scale:
    return Scale("size", "continuous", range=range, **kwargs)


def scale_alpha(range: Tuple[float, float] = (0.1, 1.0), **kwargs: Any) -> Scale:
    return Scale("alpha", "continuous", range=range, **kwargs)


def scale_shape_manual(values: Any = None, **kwargs: Any) -> Scale:
    return Scale("shape", "manual", values=values, **kwargs)


def scale_linetype_manual(values: Any = None, **kwargs: Any) -> Scale:
    return Scale("linetype", "manual", values=values, **kwargs)


# ggtree::scale_color & add_colorbar (from R/method-scale-color.R)
_DEFAULT_PALETTE = [
    "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF",
    "#E76BF3", "#00BA38", "#619CFF", "#B79F00",
    "#00B0F6", "#FF61CC",
]


def _hex_to_rgb(h: str) -> Tuple[int, int, int]:
    named = {
        "red": (255, 0, 0),
        "green": (0, 170, 0),
        "blue": (0, 0, 255),
        "black": (0, 0, 0),
        "white": (255, 255, 255),
        "grey": (128, 128, 128),
        "gray": (128, 128, 128),
        "grey50": (128, 128, 128),
        "gray50": (128, 128, 128),
        "steelblue": (70, 130, 180),
        "firebrick": (178, 34, 34),
        "darkgreen": (0, 100, 0),
        "orange": (255, 165, 0),
        "purple": (128, 0, 128),
        "yellow": (255, 215, 0),
        "pink": (255, 192, 203),
        "cyan": (0, 255, 255),
        "magenta": (255, 0, 255),
        "brown": (165, 42, 42),
        "navy": (0, 0, 128),
        "light blue": (173, 216, 230),
        "lightblue": (173, 216, 230),
    }
    s = str(h).strip().lower()
    if s in named:
        return named[s]
    if s.startswith("#") and len(s) == 7:
        return (int(s[1:3], 16), int(s[3:5], 16), int(s[5:7], 16))
    return (30, 100, 180)


def interpolate_color(val: float, vmin: float, vmax: float, low: str = "#0072B2", high: str = "#D55E00") -> str:
    if vmax <= vmin or _is_scalar_na(val):
        t = 0.5
    else:
        t = max(0.0, min(1.0, (float(val) - vmin) / (vmax - vmin)))
    r1, g1, b1 = _hex_to_rgb(low)
    r2, g2, b2 = _hex_to_rgb(high)
    r = int(round(r1 + t * (r2 - r1)))
    g = int(round(g1 + t * (g2 - g1)))
    b = int(round(b1 + t * (b2 - b1)))
    return f"#{r:02x}{g:02x}{b:02x}"


def scale_color(
    treeview: Any,
    by: str,
    low: Optional[str] = None,
    high: Optional[str] = None,
    na_color: str = "grey",
    **kwargs: Any,
) -> List[str]:
    """
    Map a tree variable 'by' to a vector of colors ordered by node ID
    (exact port of ggtree::scale_color in R/method-scale-color.R).
    """
    from .layout import fortify
    df = treeview.data if hasattr(treeview, "data") else fortify(treeview)
    df_sorted = df.sort_values("node")
    if by not in df_sorted:
        raise KeyError(f"Column '{by}' not found in tree data")
    vals = df_sorted[by].tolist()
    non_na = [v for v in vals if not _is_scalar_na(v)]
    is_num = len(non_na) > 0 and all(isinstance(v, (int, float, np.number)) and not isinstance(v, bool) for v in non_na)
    low_c = low or "#0072B2"
    high_c = high or "#D55E00"

    if is_num:
        fvals = [float(v) for v in non_na]
        vmin, vmax = min(fvals), max(fvals)
        colors = [na_color if _is_scalar_na(v) else interpolate_color(float(v), vmin, vmax, low_c, high_c) for v in vals]
        attr_dict = {"by": by, "low": low_c, "high": high_c, "min": vmin, "max": vmax}
        if hasattr(treeview, "data"):
            treeview.data.attrs["scale_color_attr"] = attr_dict
        return colors
    else:
        uniq = []
        for v in non_na:
            if v not in uniq:
                uniq.append(v)
        cmap = {u: _DEFAULT_PALETTE[i % len(_DEFAULT_PALETTE)] for i, u in enumerate(uniq)}
        return [na_color if _is_scalar_na(v) else cmap[v] for v in vals]


def add_colorbar(
    treeview: Any,
    color: Sequence[str],
    font_size: float = 4.0,
     linesize: float = 2.0,
    **kwargs: Any,
) -> Any:
    """Add a colorbar legend annotation to a ggtree plot (port of add_colorbar in R/method-scale-color.R)."""
    out = treeview.copy()
    out.theme_spec["colorbar"] = {"color": list(color), "font_size": font_size, "linesize": linesize, **kwargs}
    return out


def set_hilight_legend(p: Any, color: Sequence[str], label: Sequence[str], alpha: float = 1.0) -> Any:
    """Defunct function from R/defunct.R (geom_hilight now supports legends directly via aes(fill=...))."""
    import warnings
    warnings.warn(
        "This function is defunct. New version of geom_hilight can generate legend automatically.",
        DeprecationWarning,
        stacklevel=2,
    )
    out = p.copy()
    out.theme_spec["hilight_legend"] = {"color": list(color), "label": list(label), "alpha": float(alpha)}
    return out
