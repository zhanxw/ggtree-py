# Changelog

## 0.1.0 (2026-11-14)

Initial Python release (`ggtree-py`, imported as `ggtree`), ported from Bioconductor `ggtree` (`4.3.2`).

### Highlights
- **Tree parsers & structures (`ggtree.tree`)**: Newick (`read_tree`), BEAST NEXUS (`read_beast`), NHX (`read_nhx`), jplace (`read_jplace`), FASTA (`read_fasta`), `Phylo`, `TreeData`, `MultiPhylo`, `TreeDataList`, `Hclust`/`hclust`, and bundled sample datasets (`system_file`).
- **Tidy tree operations (`ggtree.tidytree`)**: `MRCA`, `groupOTU`, `groupClade`, `offspring`, `ancestor`, `parent`, `child`, `sibling`, `nodeid`, `nodelab`, `tree_subset`, `td_filter`, `td_mutate`, and `td_unnest`.
- **12+ tree layouts (`ggtree.layout`)**: `rectangular`, `slanted`, `circular`, `fan`, `inward_circular`, `dendrogram`, `equal_angle`, `daylight` (`unrooted`), `tidy`, `ellipse`, `roundrect`, `radial`, `ape`, and `tree_and_leaf`, plus time-scaled MRSD coordinates (`Date2decimal`, `decimal2Date`).
- **Grammar-of-graphics geoms & scales (`ggtree.aes`, `ggtree.geoms`, `ggtree.themes_scales`)**: Symbolic expressions (`aes`, `isTip`, `parent`), `geom_tree`, `geom_tiplab`, `geom_tiplab2`, `geom_nodelab`, `geom_point2`, `geom_tippoint`, `geom_nodepoint`, `geom_rootpoint`, `geom_rootedge`, `geom_aline`, `geom_hilight`, `geom_balance`, `geom_cladelab`, `geom_strip`, `geom_taxalink`, `geom_tanglelink`, `geom_range`, `geom_treescale`, `new_scale_fill`, `new_scale_color`, Viridis/Magma/Plasma/Inferno/Cividis/Brewer/gradient scales, and `theme_tree`, `theme_tree2`, `theme_dendrogram`.
- **Manipulation & composite plots (`ggtree.manipulate`, `ggtree.compose`)**: `collapse`, `expand`, `rotate`, `flip`, `scaleClade`, `viewClade`, `zoomClade`, `open_tree`, `rotate_tree`, `gheatmap`, `msaplot`, `geom_facet`/`facet_plot`, `nodepie`/`nodebar`/`inset`, `ggdensitree`, `as_tanglegram`/`ggdoubletree`, `multiplot`, `plot_list`, and `|`/`/` patchwork composition.
- **Zero-external-dependency SVG & PNG renderer (`ggtree.plot`)**: Direct vector SVG rendering and PNG export via `ggsave()`, `to_svg()`, and `to_png()` (using `rsvg-convert` when present, with pure-NumPy RGBA fallback).
- **Verified**: All 49 tests pass (`tests/test_examples.py` + `tests/test_ggtree.py`).
