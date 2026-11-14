# API Reference

All public symbols listed below are re-exported at top level (`import ggtree`) and organized into submodules.

## Core Plot & Export (`ggtree.plot`)

- `ggtree(tr, mapping=None, layout="rectangular", open_angle=0.0, mrsd=None, as_Date=False, yscale="none", yscale_mapping=None, ladderize=True, right=False, branch_length="branch.length", root_position=0.0, xlim=None, **kwargs) -> GGTree`
- `ggplot(data=None, mapping=None, **kwargs) -> GGTree`
- `ggplot_build(plot) -> BuiltPlot`
- `is_ggplot(x) -> bool`
- `ggsave(filename, plot=None, width=8, height=6, dpi=100, **kwargs) -> str`
- `last_plot() -> Optional[GGTree]`
- `GGTree.to_svg(width=800, height=600) -> str` / `GGTree.to_png(width=800, height=600) -> bytes` / `GGTree.save(...)`
- `GGTree` operators: `+` (layer/modifier), `%` / `<<` / `@` (attach metadata `DataFrame` or swap tree), `|` (horizontal compose), `/` (vertical compose)

## Tree Parsers, Structures & Generators (`ggtree.tree`)

- `Phylo`, `TreeData`, `MultiPhylo`, `TreeDataList`, `Hclust`
- `read_tree(file=None, text=None, tree_names=None, keep_multi=False)`
- `read_beast(file)` / `read_nhx(file=None, text=None)` / `read_jplace(file)` / `read_fasta(file)`
- `write_tree(phy, file=None) -> str`
- `rtree(n, rooted=True, tip_label=None, br="runif")` / `rmtree(N, n, ...)` / `stree(n, type="left", ...)` / `hclust(d, method="complete", labels=None)`
- `as_phylo(x)` / `as_treedata(x)` / `as_tibble(x) -> TblTree`
- `as_polytomy(tree, feature="node.label", fun=...)`
- `extract_clade(phy, node)` / `subtrees(tree)` / `ladderize(phy, right=True)` / `reorder_phylo(phy, order="cladewise")` / `di2multi(phy, tol=1e-8)` / `estimate_dates(tree, node_dates, mu=1.0, nsteps=1)`
- `Ntip(x)`, `Nnode(x, internal_only=True)`, `Nnode2(x)`, `getNodeNum(x)`, `is_rooted(x)`, `get_tree(x)`, `get_seqs(x)`, `system_file(*parts, package=None)`

## Tidytree Traversal & Helpers (`ggtree.tidytree`)

- `taxa2node(data, taxa)`, `nodeid(tree, label)`, `nodelab(tree, id)`, `rootnode(x)`, `getRoot(x)`, `itself(data, node)`
- `child(data, node)`, `parent(data, node)`, `sibling(data, node)`, `ancestor(data, node)`
- `offspring(data, node, tiponly=False, self_include=False, type="all")`, `offspring_bulk(data, nodes, self_include=False)`
- `MRCA(data, node1=None, node2=None, *more_nodes)`
- `groupClade(data, _node, group_name="group", overlap="overwrite")`
- `groupOTU(data, _node, group_name="group", overlap="overwrite", connect=False)`
- `tree_subset(tree, node, levels_back=5, group_node=True)`
- `td_filter(*conditions, **kwargs)`, `td_mutate(**kwargs)`, `td_unnest(cols)`

## Layouts & Coordinates (`ggtree.layout`)

- `fortify(model, data=None, layout="rectangular", ladderize=True, right=False, branch_length="branch.length", mrsd=None, as_Date=False, yscale="none", root_position=0.0)`
- `tree_layout(...)`, `get_layout(plot)`, `tip_y_alignment(x, tol=1e-6)`, `check_reverse(data)`
- `layoutEqualAngle`, `layoutDaylight`, `layoutApe`, `layoutTreeAndLeaf`, `applyLayoutDaylight`
- `Date2decimal(x)`, `decimal2Date(x)`, `scaleX_by_time_from_mrsd(df, mrsd, as_Date=False)`, `scaleY(...)`
- `get_clade_position(treeview, node)`, `get_clade_position_(data, node, reverse=False)`, `clade_geometry(...)`, `adjust_cladelabel_angle(angle, horizontal=True)`

## Aesthetics & Symbols (`ggtree.aes`)

- `aes(*args, **kwargs) -> Aes`, `Aes`, `SymExpr`, `CallableSymExpr`, `modifyList(x, val)`, `eval_expr(...)`, `eval_subset(...)`
- Predefined symbols: `isTip`, `parent`, `node`, `label`, `branch`, `branch_length`, `x`, `y`, `angle`, `group`, `rate`, `posterior`, `height`, `value`, `id`, `_id`, `place`, `type_`, `annote`, `image`, `offset`, `ta1`, `ta2`, `taxa1`, `taxa2`, `tree`, `v`, `S`, `trait`, `score`

## Geoms (`ggtree.geoms`)

- **Tree & label geoms**: `geom_tree`, `geom_tree2`, `stat_tree`, `geom_tiplab`, `geom_tiplab2`, `geom_nodelab`, `geom_nodelab2`, `geom_tippoint`, `geom_nodepoint`, `geom_rootpoint`, `geom_point2`, `geom_text2`, `geom_label2`, `geom_segment2`, `geom_aline`, `geom_rootedge`
- **Clade & phylogenetic annotation**: `geom_hilight`, `geom_highlight`, `geom_hilight_encircle`, `geom_balance`, `geom_cladelab`, `geom_cladelabel`, `geom_cladelabel2`, `geom_strip`, `geom_striplab`, `geom_taxalink`, `geom_tanglelink`, `geom_treescale`, `geom_range`, `geom_zoom_clade`, `geom_motif`, `geom_inset`
- **General / facet geoms**: `geom_point`, `geom_text`, `geom_label`, `geom_segment`, `geom_bar`, `geom_col`, `geom_tile`, `geom_raster`, `geom_line`, `geom_path`, `geom_area`, `geom_ribbon`, `geom_rect`, `geom_polygon`, `geom_curve`, `geom_vline`, `geom_hline`, `geom_abline`, `geom_boxplot`, `geom_violin`, `geom_density`, `geom_histogram`, `geom_errorbar`, `geom_errorbarh`, `geom_linerange`, `geom_pointrange`
- **Interactive aliases & helpers**: `geom_*_interactive`, `ggtree_set_interactive`, `ggtree_unset_interactive`

## Tree Manipulation (`ggtree.manipulate`)

- `collapse(x, node=None, mode="none", clade_name=None, height=None)`, `expand(tree_view, node=None)`
- `rotate(tree_view, node=None)`, `flip(tree_view, node1=None, node2=None)`
- `scaleClade(tree_view, node=None, scale=1.0, vertical_only=True)`, `viewClade(tree_view, node=None, xmax_adjust=0.0)`, `zoomClade(tree_view, node=None, xexpand=None)`
- `open_tree(p, angle=180.0)`, `rotate_tree(treeview, angle)`
- `layout_rectangular()`, `layout_circular()`, `layout_fan(angle=180.0)`, `layout_inward_circular(xlim=None)`, `layout_dendrogram()`, `layout_ggtree(layout="rectangular", angle=180.0, xlim=None)`
- `reroot(object, node)`, `rescale_tree(tree_object, branch_length="branch.length")`, `gzoom(object, focus, subtree=False, widths=(0.3, 0.7), xmax_adjust=0.0)`, `get_taxa_name(tree_view, node=None)`
- Formatting helpers: `label_pad(label, justify="right", pad=".")`, `range_format(x, digits=2)`, `scale_color_subtree` / `scale_colour_subtree`

## Composite Visualizations (`ggtree.compose`)

- `gheatmap(p, data, offset=0.0, width=1.0, low="green", high="red", color="white", colnames=True, colnames_position="bottom", colnames_angle=0.0, colnames_level=None, colnames_offset_x=0.0, colnames_offset_y=0.0, font_size=4.0, family="", hjust=0.5, legend_title="value", custom_column_labels=None)`
- `get_heatmap_column_position(treeview, by="bottom")`
- `msaplot(p, fasta, offset=0.0, width=1.0, color=None, window=None, bg_line=True, height=0.8)`
- `geom_facet(panel=None, data=None, geom=geom_tile, mapping=None, **kwargs)` / `facet_plot(...)` / `facet_data(p, panel)` / `facet_labeller(p, label)` / `facet_widths(p, widths)`
- `nodebar(data, cols, color=None, alpha=1.0, position="stack")`, `nodepie(data, cols, color=None, alpha=1.0)`, `inset(tree_view, insets, width=0.1, height=0.1, hjust=0.0, vjust=0.0, x="node", reverse_x=False, reverse_y=False)`
- `ggdensitree(data=None, mapping=None, layout="slanted", tip_order="mode", align_tips=True, jitter=0.0)`
- `Tanglegram`, `as_tanglegram(x, y=None, assoc=None, gap=0.08, optimize=False, optimize_side="right", mirror=True)`, `ggdoubletree(tr1, tr2=None, assoc=None, layout="rectangular", gap=0.08, optimize=False, optimize_side="right", mirror=True, link_color="grey50", link_alpha=0.8, ...)`
- `PlotList`, `multiplot(*plots, plotlist=None, ncol=1, widths=None, labels=None)`, `plot_list(*plots, ncol=None, nrow=None, widths=None, heights=None, tag_levels=None)`

## Themes, Scales & Coordinates (`ggtree.themes_scales`)

- **Themes & Labels**: `theme_tree`, `theme_tree2`, `theme_dendrogram`, `theme_inset`, `theme_bw`, `theme_minimal`, `theme_void`, `theme_classic`, `theme`, `labs`, `ggtitle`, `xlab`, `ylab`, `element_blank`, `element_text`, `element_line`, `element_rect`, `rel`, `unit`, `margin`, `arrow`
- **Limits & Expansion**: `xlim`, `ylim`, `xlim_tree`, `xlim_expand`, `hexpand`, `vexpand`, `ggexpand`, `revts`
- **Multi-Scale Support**: `new_scale_fill()`, `new_scale_color()`, `new_scale_colour()`
- **Scales**: `scale_fill_viridis_c`, `scale_fill_viridis_d`, `scale_color_viridis_c`, `scale_color_viridis_d`, `scale_fill_gradient`, `scale_fill_gradient2`, `scale_fill_continuous`, `scale_fill_brewer`, `scale_fill_manual`, `scale_color_gradient`, `scale_color_gradient2`, `scale_color_continuous`, `scale_color_brewer`, `scale_color_manual`, `scale_x_continuous`, `scale_y_continuous`, `scale_x_date`, `scale_x_ggtree`, `scale_x_range`, `scale_size`, `scale_alpha`, `scale_shape_manual`, `scale_linetype_manual`, `add_colorbar`, `set_hilight_legend`
- **Coords & Facets**: `coord_cartesian`, `coord_flip`, `coord_polar`, `coord_fixed`, `facet_wrap`, `facet_grid`, `guides`, `guide_legend`, `guide_colorbar`

## DataFrame & R Compatibility Helpers (`ggtree.dataframe`)

- `DataFrame`, `TblTree`, `ColumnSeries`, `Factor`
- `data_frame`, `tibble`, `factor`, `read_table`, `read_csv`
- `left_join`, `full_join`, `inner_join`, `bind_rows`, `mutate`, `filter`, `select`, `arrange`, `unnest`, `match`, `is_na`, `attr`, `set_attr`, `options`, `getOption`
- Vector/RNG/constants: `c`, `rep`, `seq`, `seq_len`, `seq_along`, `paste`, `paste0`, `rnorm`, `runif`, `rbinom`, `sample`, `set_seed`, `LETTERS`, `letters`, `NA`, `TRUE`, `FALSE`, `NULL`
- Dot-call compatibility namespaces: `read.*`, `write.*`, `system.*`, `data.*`, `set_.*`, `as_.*`, `is_.*`, `ape.*`, `treeio.*`, `tidytree.*`, `ggplot2.*`
