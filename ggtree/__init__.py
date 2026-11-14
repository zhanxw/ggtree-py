"""
ggtree: Grammar of Graphics for Phylogenetic Trees in Python (v4.3.2)
Complete Python port of the R/Bioconductor ggtree package.
"""

__version__ = "0.1.2"
__r_ggtree_version__ = "4.3.2"

from .dataframe import (
    DataFrame,
    TblTree,
    ColumnSeries,
    Factor,
    options,
    getOption,
    attr,
    set_attr,
    data_frame,
    tibble,
    c,
    rep,
    seq,
    seq_len,
    seq_along,
    paste,
    paste0,
    factor,
    rnorm,
    runif,
    rbinom,
    sample,
    set_seed,
    LETTERS,
    letters,
    NA,
    TRUE,
    FALSE,
    NULL,
    is_na,
    match,
    mutate,
    filter,
    select,
    arrange,
    left_join,
    full_join,
    inner_join,
    bind_rows,
    unnest,
    read_table,
    read_csv,
)

from .tree import (
    Phylo,
    TreeData,
    MultiPhylo,
    TreeDataList,
    Hclust,
    hclust,
    read_tree,
    read_beast,
    read_nhx,
    read_jplace,
    read_fasta,
    write_tree,
    rtree,
    rmtree,
    stree,
    estimate_dates,
    reorder_phylo,
    ladderize,
    di2multi,
    as_polytomy,
    extract_clade,
    Ntip,
    Nnode,
    Nnode2,
    getNodeNum,
    is_rooted,
    get_tree,
    as_phylo,
    as_treedata,
    as_tibble,
    system_file,
    get_seqs,
    subtrees,
)

from .tidytree import (
    taxa2node,
    nodeid,
    nodelab,
    rootnode,
    getRoot,
    isTip as _tidytree_isTip,
    itself,
    child,
    parent as _tidytree_parent,
    sibling,
    ancestor,
    offspring,
    offspring_bulk,
    MRCA,
    groupOTU,
    groupClade,
    tree_subset,
    td_filter,
    td_mutate,
    td_unnest,
)

from .layout import (
    Date2decimal,
    decimal2Date,
    getXcoord2,
    getXcoord,
    getXcoord_no_length,
    getXcoord_no_length_slanted,
    getYcoord_no_length_slanted,
    getYcoord,
    getYcoord_tidy,
    re_assign_ycoord_df,
    reassign_y_from_node_to_root,
    calculate_branch_mid,
    calculate_angle,
    add_angle_slanted,
    set_branch_length_cladogram,
    getNodes_by_postorder,
    getNodesBreadthFirst_df,
    layoutEqualAngle,
    applyLayoutDaylight,
    layoutDaylight,
    layoutApe,
    layoutTreeAndLeaf,
    scaleX_by_time_from_mrsd,
    scaleY,
    adjust_hclust_tip_edge_len,
    tree_layout,
    fortify,
    get_layout,
    tip_y_alignment,
    check_reverse,
    get_clade_position_,
    get_clade_position,
    clade_geometry,
    adjust_cladelabel_angle,
)

from .aes import (
    SymExpr,
    CallableSymExpr,
    Aes,
    aes,
    modifyList,
    eval_expr,
    eval_subset,
    node,
    label,
    branch,
    branch_length,
    x,
    y,
    angle,
    group,
    rate,
    posterior,
    height,
    value,
    id,
    _id,
    place,
    type_,
    annote,
    image,
    offset,
    ta1,
    ta2,
    taxa1,
    taxa2,
    tree,
    v,
    S,
    trait,
    score,
)

isTip = CallableSymExpr("isTip", _tidytree_isTip)
parent = CallableSymExpr("parent", _tidytree_parent)

from .themes_scales import (
    PlotModifier,
    Theme,
    Scale,
    Coord,
    rel,
    unit,
    margin,
    arrow,
    element_blank,
    element_text,
    element_line,
    element_rect,
    guide_legend,
    guide_colorbar,
    guides,
    theme,
    theme_tree,
    theme_tree2,
    theme_dendrogram,
    theme_inset,
    theme_bw,
    theme_minimal,
    theme_void,
    theme_classic,
    labs,
    ggtitle,
    xlab,
    ylab,
    xlim,
    ylim,
    xlim_expand,
    xlim_tree,
    hexpand,
    vexpand,
    ggexpand,
    revts,
    scale_x_range,
    coord_cartesian,
    coord_flip,
    coord_polar,
    coord_fixed,
    facet_wrap,
    facet_grid,
    new_scale_fill,
    new_scale_color,
    new_scale_colour,
    scale_color_manual,
    scale_colour_manual,
    scale_fill_manual,
    scale_color_continuous,
    scale_colour_continuous,
    scale_fill_continuous,
    scale_color_gradient,
    scale_colour_gradient,
    scale_fill_gradient,
    scale_color_gradient2,
    scale_fill_gradient2,
    scale_color_viridis_c,
    scale_colour_viridis_c,
    scale_fill_viridis_c,
    scale_color_viridis_d,
    scale_colour_viridis_d,
    scale_fill_viridis_d,
    scale_color_brewer,
    scale_fill_brewer,
    scale_x_continuous,
    scale_y_continuous,
    scale_x_date,
    scale_x_ggtree,
    scale_size,
    scale_size_continuous,
    scale_alpha,
    scale_shape_manual,
    scale_linetype_manual,
    scale_color,
    add_colorbar,
    set_hilight_legend,
)

from .geoms import (
    Layer,
    LayerList,
    geom_tree,
    geom_tree2,
    stat_tree,
    geom_tiplab,
    geom_tiplab2,
    geom_nodelab,
    geom_nodelab2,
    geom_point2,
    geom_point,
    geom_tippoint,
    geom_nodepoint,
    geom_rootpoint,
    geom_text2,
    geom_text,
    geom_label2,
    geom_label,
    geom_segment2,
    geom_segment,
    geom_aline,
    geom_hilight,
    geom_highlight,
    geom_hilight_encircle,
    geom_cladelab,
    geom_cladelabel,
    geom_cladelabel2,
    geom_strip,
    geom_striplab,
    geom_taxalink,
    geom_tanglelink,
    geom_treescale,
    geom_balance,
    geom_range,
    geom_rootedge,
    geom_zoom_clade,
    geom_motif,
    geom_inset,
    geom_bar,
    geom_col,
    geom_tile,
    geom_raster,
    geom_line,
    geom_path,
    geom_area,
    geom_ribbon,
    geom_rect,
    geom_polygon,
    geom_curve,
    geom_vline,
    geom_hline,
    geom_abline,
    geom_boxplot,
    geom_violin,
    geom_density,
    geom_histogram,
    geom_errorbar,
    geom_errorbarh,
    geom_linerange,
    geom_pointrange,
    geom_tree_interactive,
    geom_point_interactive,
    geom_tiplab_interactive,
    geom_tiplab2_interactive,
    geom_nodelab_interactive,
    geom_nodelab2_interactive,
    geom_cladelab_interactive,
    geom_cladelabel_interactive,
    geom_hilight_interactive,
    geom_strip_interactive,
    geom_taxalink_interactive,
    geom_label2_interactive,
    geom_point2_interactive,
    geom_segment2_interactive,
    geom_text2_interactive,
    GeomHilightEncircle,
    GeomInteractiveCurvelink,
    GeomInteractiveHilightEncircle,
    GeomInteractiveHilightRect,
    GeomInteractiveLabelGGtree,
    GeomInteractivePointGGtree,
    GeomInteractiveSegmentGGtree,
    GeomInteractiveTextGGtree,
    StatBalance,
    ggtree_set_interactive,
    ggtree_unset_interactive,
)

from .manipulate import (
    get_taxa_name,
    viewClade,
    zoomClade,
    collapse,
    expand,
    rotate,
    flip,
    scaleClade,
    open_tree,
    rotate_tree,
    layout_ggtree,
    layout_rectangular,
    layout_circular,
    layout_fan,
    layout_inward_circular,
    layout_dendrogram,
    label_pad,
    range_format,
    reroot,
    rescale_tree,
    gzoom,
    identify,
    scale_color_subtree,
    scale_colour_subtree,
)

from .compose import (
    gheatmap,
    get_heatmap_column_position,
    msaplot,
    geom_facet,
    facet_plot,
    facet_data,
    facet_labeller,
    facet_widths,
    nodebar,
    nodepie,
    inset,
    ggdensitree,
    Tanglegram,
    as_tanglegram,
    as_tanglegram_cophylo,
    as_tanglegram_default,
    as_tanglegram_tanglegram,
    ggdoubletree,
    PlotList,
    multiplot,
    plot_list,
)

from .plot import (
    GGTree,
    BuiltPlot,
    BuiltLayout,
    ggplot,
    ggtree,
    ggplot_build,
    is_ggplot,
    ggsave,
    last_plot,
)


# ---------------------------------------------------------------------------
# R-namespace dot-Call Compatibility Proxies
# Allows R-style expressions in Python:
#   read.tree(...), read.beast(...), system.file(...), data.frame(...),
#   set.seed(...), as.polytomy(...), ape.rtree(...), ggplot2.ggplot_build(...)
# ---------------------------------------------------------------------------

class _ReadNamespace:
    tree = staticmethod(read_tree)
    beast = staticmethod(read_beast)
    nhx = staticmethod(read_nhx)
    jplace = staticmethod(read_jplace)
    fasta = staticmethod(read_fasta)
    table = staticmethod(read_table)
    csv = staticmethod(read_csv)


class _WriteNamespace:
    tree = staticmethod(write_tree)


class _SystemNamespace:
    file = staticmethod(system_file)


class _DataNamespace:
    frame = staticmethod(data_frame)


class _SetNamespace:
    seed = staticmethod(set_seed)


class _AsNamespace:
    polytomy = staticmethod(as_polytomy)
    phylo = staticmethod(as_phylo)
    treedata = staticmethod(as_treedata)
    tibble = staticmethod(as_tibble)
    tanglegram = staticmethod(as_tanglegram)


class _IsNamespace:
    na = staticmethod(is_na)
    ggplot = staticmethod(is_ggplot)
    rooted = staticmethod(is_rooted)


class _ApeNamespace:
    rtree = staticmethod(rtree)
    rmtree = staticmethod(rmtree)
    stree = staticmethod(stree)
    read_tree = staticmethod(read_tree)
    write_tree = staticmethod(write_tree)
    read = _ReadNamespace()
    write = _WriteNamespace()
    Ntip = staticmethod(Ntip)
    Nnode = staticmethod(Nnode)
    di2multi = staticmethod(di2multi)
    ladderize = staticmethod(ladderize)
    extract_clade = staticmethod(extract_clade)
    estimate_dates = staticmethod(estimate_dates)
    is_rooted = staticmethod(is_rooted)
    subtrees = staticmethod(subtrees)


class _TreeioNamespace:
    read_tree = staticmethod(read_tree)
    read_beast = staticmethod(read_beast)
    read_nhx = staticmethod(read_nhx)
    read_jplace = staticmethod(read_jplace)
    read_fasta = staticmethod(read_fasta)
    read = _ReadNamespace()
    get_tree = staticmethod(get_tree)
    get_seqs = staticmethod(get_seqs)
    Nnode = staticmethod(Nnode)
    Nnode2 = staticmethod(Nnode2)
    Ntip = staticmethod(Ntip)
    as_phylo = staticmethod(as_phylo)
    as_treedata = staticmethod(as_treedata)
    rescale_tree = staticmethod(rescale_tree)


class _TidytreeNamespace:
    nodeid = staticmethod(nodeid)
    nodelab = staticmethod(nodelab)
    rootnode = staticmethod(rootnode)
    isTip = staticmethod(_tidytree_isTip)
    child = staticmethod(child)
    parent = staticmethod(_tidytree_parent)
    sibling = staticmethod(sibling)
    ancestor = staticmethod(ancestor)
    offspring = staticmethod(offspring)
    MRCA = staticmethod(MRCA)
    groupOTU = staticmethod(groupOTU)
    groupClade = staticmethod(groupClade)
    tree_subset = staticmethod(tree_subset)
    as_tibble = staticmethod(as_tibble)
    as_treedata = staticmethod(as_treedata)
    as_phylo = staticmethod(as_phylo)


class _Ggplot2Namespace:
    ggplot = staticmethod(ggplot)
    aes = staticmethod(aes)
    ggplot_build = staticmethod(ggplot_build)
    is_ggplot = staticmethod(is_ggplot)
    ggsave = staticmethod(ggsave)
    rel = staticmethod(rel)
    unit = staticmethod(unit)
    margin = staticmethod(margin)
    arrow = staticmethod(arrow)
    element_blank = staticmethod(element_blank)
    element_text = staticmethod(element_text)
    element_line = staticmethod(element_line)
    element_rect = staticmethod(element_rect)
    theme = staticmethod(theme)
    theme_bw = staticmethod(theme_bw)
    theme_minimal = staticmethod(theme_minimal)
    theme_void = staticmethod(theme_void)
    theme_classic = staticmethod(theme_classic)
    labs = staticmethod(labs)
    ggtitle = staticmethod(ggtitle)
    xlab = staticmethod(xlab)
    ylab = staticmethod(ylab)
    xlim = staticmethod(xlim)
    ylim = staticmethod(ylim)
    coord_cartesian = staticmethod(coord_cartesian)
    coord_flip = staticmethod(coord_flip)
    coord_polar = staticmethod(coord_polar)
    coord_fixed = staticmethod(coord_fixed)
    facet_wrap = staticmethod(facet_wrap)
    facet_grid = staticmethod(facet_grid)
    fortify = staticmethod(fortify)


read = _ReadNamespace()
write = _WriteNamespace()
system = _SystemNamespace()
data = _DataNamespace()
set_ = _SetNamespace()
as_ = _AsNamespace()
is_ = _IsNamespace()
ape = _ApeNamespace()
treeio = _TreeioNamespace()
tidytree = _TidytreeNamespace()
ggplot2 = _Ggplot2Namespace()
