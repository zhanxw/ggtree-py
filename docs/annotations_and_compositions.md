# Annotations, Heatmaps, Facets & Multi-Tree Compositions

> **Navigation:** [Docs Home](index.md) · [Getting Started](getting_started.md) · [Layouts & Manipulation](layouts_and_manipulation.md) · [Annotations & Compositions](annotations_and_compositions.md) · [API Reference](api_reference.md) · [R-to-Python Migration](r_to_python_migration.md) · [GitHub](https://github.com/zhanxw/ggtree.py) · [PyPI](https://pypi.org/project/ggtree-py/)

## 1. Tree Annotation Layers (`ggtree.geoms`)

| Geom | Purpose |
| --- | --- |
| `geom_tree(...)`, `geom_tree2(...)` | Render tree edges (horizontal/vertical or slanted/unrooted segments) |
| `geom_tiplab(align=False, linetype="dotted", offset=0.0, hjust=0.0, ...)` | Tip labels with optional dotted/dashed alignment leader lines |
| `geom_tiplab2(...)` | Circular-oriented tip labels (auto-flips text on left semicircle) |
| `geom_nodelab(nudge_x=0.0, nudge_y=0.0, geom="text", hjust=0.5, node="internal")` | Internal (or `"all"`, `"external"`) node labels |
| `geom_nodelab2(...)` | Circular-oriented node labels |
| `geom_tippoint(...)`, `geom_nodepoint(...)`, `geom_rootpoint(...)` | Filtered point layers (`isTip`, `!isTip`, `node == parent`) |
| `geom_point2(...)`, `geom_text2(...)`, `geom_label2(...)`, `geom_segment2(...)` | Subset-aware (`aes(subset=...)`) general layers |
| `geom_rootedge(rootedge=None)` | Draw root branch stem on left side of root |
| `geom_aline(...)` | Leader lines connecting tips to aligned annotations |
| `geom_hilight(node=..., fill="steelblue", alpha=0.5, extend=0.0, to_bottom=False, type="auto")` | Highlight rectangular or polar clade wedges |
| `geom_balance(node=..., fill="steelblue", ...)` | Highlight two sibling subclades with distinct fill bands |
| `geom_cladelab(node=..., label=..., offset=0.0, offset_text=0.0, align=False, parse=False, ...)` | Bar + text clade annotation (`geom_cladelabel` / `geom_cladelabel2` aliases) |
| `geom_strip(taxa1=..., taxa2=..., label=..., offset=0.0, angle=0.0, align=True, ...)` | Span bar + label between two arbitrary taxa (`geom_striplab` alias; supports `angle="auto"`) |
| `geom_taxalink(taxa1=..., taxa2=..., curvature=0.5, outward="auto", ...)` | Curved link arc between pairs of taxa |
| `geom_range(range="length_0.95_HPD", center="auto", color="red", size=2.0, alpha=0.6)` | Horizontal HPD / confidence interval bars (supports both `_0.95_HPD` and `_95%_HPD`) |
| `geom_treescale(x=None, y=None, width=None, ...)` | Evolutionary distance scale bar with numeric width label |
| `geom_motif(on=None, label=None, align="left")` | Gene/domain motif track geometry |

`parse=True` supports basic R plotmath strings such as `italic('Clade IL')`, `bold('Label')`, and `bolditalic('Taxon')`.

---

## 2. Aligned Heatmaps & Sequences (`gheatmap`, `msaplot`, `new_scale_fill`)

### `gheatmap()`

Appends a matrix/DataFrame heatmap aligned row-for-row with tree tips (works in both rectangular and circular/polar layouts):

```python
p_hm = gg.gheatmap(
    p_tree,
    genotype_df,
    offset=5.0,
    width=0.5,
    font_size=3.0,
    colnames_angle=-45,
    hjust=0.0,
)
```

- Column positions can be inspected with `get_heatmap_column_position(p_hm, by="bottom")`.
- Missing values (`NA`, `None`, `""`) are handled cleanly.

### Multiple Fill/Color Scales (`new_scale_fill()`, `new_scale_color()`)

Track independent color/fill legends across stacked `gheatmap()` calls:

```python
p1 = gg.gheatmap(base, df_cat, offset=0.8, width=0.2) + gg.scale_fill_viridis_d(option="D") + gg.new_scale_fill()
p2 = gg.gheatmap(p1, df_num, offset=12.8, width=0.3) + gg.scale_fill_viridis_c(option="A")
```

Supported continuous & discrete palettes:
- `scale_fill_viridis_c` / `scale_fill_viridis_d` / `scale_color_viridis_c` / `scale_color_viridis_d` (options `"D"`, `"A"`, `"B"`, `"C"`, `"E"`)
- `scale_fill_gradient`, `scale_fill_gradient2`, `scale_fill_continuous`, `scale_fill_brewer`, `scale_fill_manual` (and matching `scale_color_*` / `scale_colour_*`)

### `msaplot()`

Aligns a FASTA sequence alignment next to a tree:

```python
p_msa = gg.msaplot(
    gg.ggtree(pa_tree) + gg.geom_tiplab(size=2.5),
    fasta=gg.system_file("extdata/pa.fas", package="ggtree"),
    window=(1, 60),
    offset=0.1,
    width=1.2,
)
```

---

## 3. Faceted Data Panels (`geom_facet` / `facet_plot`)

Attach arbitrary ggplot-style panels synchronized by tip `y` coordinate alongside the `"Tree"` panel:

```python
p_facet = (
    gg.ggtree(tr)
    + gg.geom_tiplab()
    + gg.geom_facet(panel="Trait", data=trait_df, geom=gg.geom_point, mapping=gg.aes(x="value"))
    + gg.xlim_tree(6)
)
p_facet = gg.facet_widths(p_facet, widths=[1, 2])
p_facet = gg.facet_labeller(p_facet, {"Tree": "Phylogeny", "Trait": "Continuous Trait"})
```

- Retrieve panel data via `facet_data(p_facet, "Trait")`.
- Supported panel geoms include `geom_point`, `geom_tile`, `geom_bar`, `geom_col`, `geom_boxplot`, `geom_violin`, `geom_density`, `geom_histogram`, `geom_segment`, `geom_line`, `geom_errorbar`, and `geom_text`.

---

## 4. Node Pie/Bar Insets (`nodepie`, `nodebar`, `inset`)

Overlay pie charts or bar charts at internal/tip nodes:

```python
anc_df = gg.DataFrame({
    "node": [11, 12, 13],
    "A": [0.7, 0.2, 0.1],
    "B": [0.2, 0.5, 0.3],
    "C": [0.1, 0.3, 0.6],
})
pies = gg.nodepie(anc_df, cols=[1, 2, 3])  # or cols=["A", "B", "C"]
p_pie = gg.inset(gg.ggtree(tr) + gg.geom_tiplab(), pies, width=0.12, height=0.12)
```

`nodebar(data, cols=..., color=None, alpha=1.0, position="stack")` provides grouped/stacked bar insets.

---

## 5. Multi-Tree Plots & Tanglegrams

### `ggdensitree()`

Overlays an ensemble of `MultiPhylo` / `TreeDataList` trees (`tip_order="mode"`, `align_tips=True`, `alpha`, `colour`):

```python
trees = gg.read_tree(gg.system_file("extdata", "ggdensitree_example.tree", package="ggtree"))
p_dens = gg.ggdensitree(trees, alpha=0.3, colour="steelblue") + gg.geom_tiplab(size=2.5) + gg.hexpand(0.3)
```

### `as_tanglegram()` & `ggdoubletree()`

Align two trees face-to-face with optional iterative crossing minimization (`optimize=True`, `optimize_side="left" | "right" | "both"`) and cubic Bezier association curves:

```python
p_tg = gg.ggdoubletree(
    gg.ggtree(t1) + gg.geom_tiplab(),
    gg.ggtree(t2) + gg.geom_tiplab(),
    assoc=assoc_df,
    optimize=True,
    optimize_side="both",
    link_color="steelblue",
)
```

### Multi-Panel Grids (`multiplot`, `plot_list`, `|`, `/`)

- `multiplot(p1, p2, ncol=2)` or `plot_list(*plots, ncol=..., nrow=..., widths=..., heights=..., tag_levels="A")`
- Nested patchwork operators: `(p1 | p2) / p3` builds a 2-column top row above a full-width bottom plot.
