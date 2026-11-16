# Visual Gallery (`ggtree-py`)

> **Navigation:** [Docs Home](index.md) · [Online Book](https://zhanxw.github.io/ggtree-py-book/) · [Visual Gallery](gallery.md) · [Getting Started](getting_started.md) · [Layouts & Manipulation](layouts_and_manipulation.md) · [Annotations & Compositions](annotations_and_compositions.md) · [API Reference](api_reference.md) · [R-to-Python Migration](r_to_python_migration.md) · [GitHub](https://github.com/zhanxw/ggtree-py) · [PyPI](https://pypi.org/project/ggtree-py/)

> 📖 **Interactive Online Book Edition:** Browse all **111 reproducible `treedata-book` figures** with live search, Part/Chapter filters, Master Index Table, and collapsible Python code in the **[Online Book Gallery](https://zhanxw.github.io/ggtree-py-book/gallery.html)** or start from **[Book Home](https://zhanxw.github.io/ggtree-py-book/)**.

All 28 rendered outputs below were generated directly by `ggtree-py` using [`docs/generate_doc_figures.py`](https://github.com/zhanxw/ggtree-py/blob/main/docs/generate_doc_figures.py) (`dpi=120`).

---

## 1. Core Tree Inputs, Metadata & Themes

See **[Getting Started & Tree Data I/O](getting_started.md)** for full code examples.

### Basic Tree Structures (`read_tree`, `branch_length="none"`, `hclust`)
*Phylogram, equal-branch cladogram, and average-linkage hierarchical clustering dendrogram rendered side by side with `|`.*

![Basic tree structures](figures/gs_01_basic_trees.png)

### Attaching Metadata via `%` (`p % df`)
*Left-joining external sample scores and groups onto tree nodes maps point size (`aes(size="score")`) and discrete Brewer colors (`scale_color_brewer`).*

![Metadata attachment via %](figures/gs_02_metadata_operator.png)

### Table Filtering & Label Formatting (`td_filter`, `td_mutate`, `label_pad`)
*Selecting tips vs internal nodes dynamically inside layer `data=` callables and right-justifying labels with `label_pad()`.*

![td_filter and td_mutate](figures/gs_03_td_filter_mutate.png)

### Symbolic Aesthetics & Branch Subsetting (`~isTip`, `node`, `branch`)
*Unquoted R-compatible expressions highlight internal node IDs and annotate non-zero branch lengths along horizontal midpoints.*

![Symbolic aesthetics and subsetting](figures/gs_04_aes_symbolic_subset.png)

### Built-In Plot Themes (`theme_tree`, `theme_tree2`, `theme_bw`, `theme_minimal`)
*Switching cleanly between blank phylogenetic canvases, evolutionary x-axis rulers, framed grids, and minimal borderless grids.*

![Themes comparison](figures/gs_05_themes_comparison.png)

---

## 2. Phylogenetic Layouts & Timescales

See **[Layouts & Tree Manipulation](layouts_and_manipulation.md)** for details.

### All 12+ Phylogenetic Layouts (`ggtree(..., layout=...)`)
*Twelve built-in layouts: `rectangular`, `roundrect`, `slanted`, `ellipse`, `circular`, `fan`, `inward_circular`, `radial`, `dendrogram`, `equal_angle`, `daylight`, and `ape`.*

![All 12 phylogenetic layouts](figures/01_all_layouts_grid.png)

### Circular Opening & Rotation (`open_tree`, `rotate_tree`)
*Progressive transformation from closed `circular` (`360°`) to `90°` fan gap via `open_tree()`, followed by `45°` counter-clockwise rotation via `rotate_tree()`.*

![Circular, open_tree, rotate_tree](figures/02_circular_open_rot.png)

### Dendrogram Layouts (`layout="dendrogram"`, `hclust`, `theme_dendrogram`)
*Vertical top-to-bottom dendrograms for clade-grouped phylogenies and distance-matrix hierarchical clusters.*

![Dendrogram layouts](figures/03_dendrogram.png)

### Time-Scaled FluA H3 Trees (`mrsd` & `revts`)
*Calibrating branch coordinates to calendar years (`mrsd="2013-01-01"`) vs reversed time-before-present (`revts() + scale_x_continuous(labels=abs)`).*

![Time-scaled mrsd and revts](figures/lm_04_mrsd_and_revts.png)

---

## 3. Clade Grouping & Tree Manipulation

See **[Layouts & Tree Manipulation](layouts_and_manipulation.md#3-clade-selection-grouping--tree-manipulation)** for API notes.

### Grouping Clades & Taxa (`groupClade`, `groupOTU`)
*Coloring monophyletic clades (`groupClade`) vs polyphyletic tip groups and connecting ancestor lineages (`groupOTU`).*

![groupClade and groupOTU](figures/lm_05_group_clade_otu.png)

### Collapsing & Expanding Clades (`collapse`, `expand`)
*Comparing `mode="mixed"`, `"max"`, and `"min"` triangle styles alongside live clade restoration via `expand()`.*

![Collapsing and expanding clades](figures/10_collapsed_and_expanded.png)

### Rotating, Flipping & Rescaling Clades (`rotate`, `flip`, `scaleClade`)
*Reordering child branches (`rotate`), swapping two sibling clades vertically (`flip`), and compressing vertical clade spacing (`scaleClade`).*

![Rotate, flip, and scaleClade](figures/lm_07_rotate_flip_scale.png)

### Focusing & Zooming Clades (`viewClade`, `gzoom`)
*Zooming viewport bounds directly onto a clade (`viewClade`) or displaying synchronized overview + focal panes (`gzoom`).*

![viewClade and gzoom](figures/lm_08_viewclade_gzoom.png)

### Collapsing Weak Branches to Polytomies (`as_polytomy`)
*Converting internal edges below bootstrap support `< 70` in `RAxML_bipartitions.H3` into clean polytomies.*

![Collapsing low-support nodes to polytomies](figures/lm_09_as_polytomy.png)

---

## 4. Annotations, Heatmaps & Sequence Alignments

See **[Annotations, Heatmaps, Facets & Multi-Tree Compositions](annotations_and_compositions.md)** for full recipes.

### Clade Highlights, Labels, Strips & Links (`geom_hilight`, `geom_cladelab`, `geom_strip`, `geom_taxalink`, `geom_treescale`)
*Combined annotation layers featuring `parse=True` plotmath formatting (`italic(...)`, `bold(...)`), root stem, and evolutionary distance bar.*

![Richly annotated tree](figures/04_annotated_tree.png)

### Clade Highlights vs Sibling Balance (`geom_hilight`, `geom_balance`)
*Comparing single-clade rectangular highlights with two-color sibling subclade balance boxes beneath `MRCA("A", "H")`.*

![Clade highlights and balance](figures/ac_02_hilight_and_balance.png)

### Span Strips & Curved Taxa Links (`geom_strip`, `geom_taxalink`)
*Annotation strips between non-monophyletic taxa bounds and Bezier curves connecting pairs across rectangular and circular trees.*

![Taxa links and strips](figures/ac_03_taxalink_and_strip.png)

### BEAST Time Tree + Uncertainty Bars + Genotype Matrix (`read_beast`, `geom_range`, `gheatmap`)
*95% HPD posterior uncertainty bars overlaid on FluA H3 nodes beside an 8-segment viral genotype heatmap (`Genotype.txt`).*

![BEAST phylogeny with genotype heatmap](figures/05_beast_gheatmap.png)

### Circular Tree with Multiple Fill Scales (`new_scale_fill`, `scale_fill_viridis_d`, `scale_fill_viridis_c`)
*Concentric categorical (`option="D"`) and continuous (`option="A"`) heatmap rings around a polar tree with independent legends.*

![Circular tree with stacked heatmaps](figures/06_circular_stacked_heatmap.png)

### Phylogeny + Multiple Sequence Alignment (`msaplot`)
*Aligned amino-acid/nucleotide FASTA sequences (`pa.nwk` + `pa.fas`, sites `1..60`) synchronized row-by-row with tip order.*

![Multiple sequence alignment plot](figures/07_msaplot.png)

---

## 5. Facets, Node Insets & Multi-Tree Compositions

See **[Annotations & Compositions](annotations_and_compositions.md#3-faceted-data-panels-geom_facet--facet_plot)** for usage details.

### Single Trait Facet Panel (`geom_facet`, `facet_widths`, `facet_labeller`)
*Synchronizing an external continuous trait scatter panel row-by-row with phylogenetic tips.*

![Single facet panel](figures/08_facet_plot.png)

### Multi-Panel Facets (`geom_col`, `geom_segment`, `geom_point`)
*Three aligned panels displaying phylogeny, habitat abundance bars (`geom_col`), and trait lollipops (`geom_segment` + `geom_point`) with a unified legend.*

![Multi-panel facets with geom_col and dotplot](figures/ac_07_multifacet_panels.png)

### Ancestral State Pie Insets (`nodepie`, `inset`)
*Overlaying discrete ancestral state proportion pie charts directly onto internal tree nodes.*

![Node pie insets](figures/09_inset_pies.png)

### Stacked vs Dodged Node Bar Insets (`nodebar`, `inset`)
*Comparing stacked proportion bars (`position="stack"`) and grouped multi-series bars (`position="dodge"`) anchored at internal nodes.*

![Stacked vs dodged nodebar insets](figures/ac_09_inset_bars.png)

### Multi-Tree Ensemble Overlay (`ggdensitree`)
*Visualizing topological uncertainty across a `MultiPhylo` posterior sample with aligned tip consensus order.*

![Densitree overlay](figures/11_ggdensitree.png)

### Tanglegram Layout & Crossing Optimization (`ggdoubletree`)
*Comparing raw face-to-face tree tip crossings (`optimize=False`) against bidirectional untangling (`optimize=True, optimize_side="both"`).*

![Tanglegram before and after optimization](figures/ac_11_tanglegram_before_after.png)

### Multi-Panel Grid Layout (`multiplot`)
*Arranging rectangular, slanted, and circular views in an indexed multi-column grid (`labels=["A", "B", "C"]`).*

![Multiplot grid](figures/13_multiplot.png)

### Nested Patchwork Composition (`|` and `/`)
*Declarative composite figures built via `(p_rect | p_circ) / p_bot`, combining top-row panels over a full-width heatmap view.*

![Patchwork composition](figures/14_patchwork.png)
