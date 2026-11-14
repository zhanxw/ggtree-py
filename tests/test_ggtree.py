"""
Comprehensive test suite porting all 18 R testthat files in `tests/testthat/test-*.R` to Python unittest.
"""

import math
import tempfile
import unittest
import warnings
from pathlib import Path
import numpy as np

from ggtree import (
    ggtree,
    ggplot,
    ggplot_build,
    is_ggplot,
    ggsave,
    fortify,
    read_tree,
    read_beast,
    read_nhx,
    read_jplace,
    read_fasta,
    rtree,
    rmtree,
    stree,
    hclust,
    as_phylo,
    as_treedata,
    as_tibble,
    as_polytomy,
    di2multi,
    extract_clade,
    estimate_dates,
    system_file,
    DataFrame,
    data_frame,
    tibble,
    c,
    rep,
    seq,
    paste0,
    rnorm,
    runif,
    sample,
    set_seed,
    LETTERS,
    letters,
    NA,
    TRUE,
    FALSE,
    is_na,
    options,
    attr,
    aes,
    isTip,
    node,
    parent,
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
    place,
    score,
    MRCA,
    nodeid,
    nodelab,
    rootnode,
    child,
    ancestor,
    offspring,
    sibling,
    groupOTU,
    groupClade,
    tree_subset,
    td_filter,
    td_mutate,
    td_unnest,
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
    layout_rectangular,
    layout_circular,
    layout_fan,
    layout_inward_circular,
    layout_dendrogram,
    label_pad,
    range_format,
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
    geom_vline,
    geom_hline,
    theme,
    theme_tree,
    theme_tree2,
    theme_dendrogram,
    theme_inset,
    theme_bw,
    theme_minimal,
    theme_void,
    rel,
    unit,
    margin,
    arrow,
    element_blank,
    element_text,
    element_line,
    element_rect,
    guide_legend,
    guides,
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
    scale_color_manual,
    scale_colour_manual,
    scale_fill_manual,
    scale_color_continuous,
    scale_fill_viridis_d,
    scale_x_ggtree,
    scale_color,
    add_colorbar,
    new_scale_fill,
    new_scale_color,
    coord_cartesian,
    coord_polar,
    coord_flip,
    facet_wrap,
    facet_grid,
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
    as_tanglegram,
    ggdoubletree,
    multiplot,
    plot_list,
    Date2decimal,
    decimal2Date,
    get_clade_position,
    clade_geometry,
    Aes,
    SymExpr,
    TreeDataList,
    left_join,
    full_join,
    factor,
    get_layout,
    tip_y_alignment,
    new_scale_colour,
    scale_colour_gradient,
    scale_colour_viridis_c,
    scale_fill_viridis_c,
    _id,
    type_,
    ape,
    ggplot2,
    treeio,
    tidytree,
    read,
    system,
)


class TestGGTreeSuite(unittest.TestCase):
    """Port of all 18 testthat files in tests/testthat/test-*.R."""

    def test_group_clade_and_group_otu(self):
        """Port of tests/testthat/test-group.R."""
        # groupClade works for collapses node
        nwk = system_file("extdata", "sample.nwk", package="treeio")
        tree = read_tree(nwk)
        p = ggtree(tree) + geom_tiplab()
        p = collapse(p, node=21)
        p = groupClade(p, 21)
        self.assertIn("group", p.data.columns)
        self.assertGreaterEqual(len(set(p.data["group"])), 2)

        # groupOTU and groupClade on phylo, treedata, tbl_tree, and ggtree
        tr = read_tree(text="((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
        tr_g = groupOTU(tr, ["D", "E"], group_name="grp")
        df_g = fortify(tr_g)
        self.assertIn("grp", df_g.columns)

        mrca_de = MRCA(tr, "D", "E")
        tr_c = groupClade(tr, mrca_de, group_name="clade_grp")
        df_c = fortify(tr_c)
        self.assertIn("clade_grp", df_c.columns)
        off = offspring(df_c, mrca_de, self_include=True)
        for nid in off["node"]:
            row = df_c.filter(lambda r: r["node"] == nid)
            self.assertEqual(str(row["clade_grp"][0]), "1")

    def test_xlim_expand(self):
        """Port of tests/testthat/test-xlim_expand.R."""
        set_seed(2019 - 5 - 9)
        x_tr = rtree(30)
        p = ggtree(x_tr) + geom_tiplab()
        d = DataFrame({"label": x_tr.tip.label, "value": rnorm(30)})
        p2 = (
            p
            + geom_facet(panel="Dot", data=d, geom=geom_point, mapping=aes(x=value))
            + xlim_tree(6)
            + xlim_expand(c(-10, 10), "Dot")
        )
        self.assertEqual(p2.panel_xlims.get("Tree"), (0.0, 6.0))
        self.assertEqual(p2.panel_xlims.get("Dot"), (-10.0, 10.0))

    def test_tidy_layout_exact_reference(self):
        """Port of tests/testthat/test-tidy-layout.R."""
        tr = read_tree(
            text="((t1:0.1,t2:0.2,(t3:0.05,t4:0.15):0.1):0.3,((t5:0.2,t6:0.1):0.05,t7:0.4):0.2,(t8:0.15,t9:0.25,t10:0.1):0.35);"
        )
        tidy_df = fortify(tr, layout="tidy")
        self.assertEqual(list(tidy_df.columns)[:9], ["parent", "node", "branch.length", "label", "isTip", "x", "y", "branch", "angle"])
        self.assertEqual(len(tidy_df), 16)
        self.assertEqual(sum(bool(v) for v in tidy_df["isTip"]), 10)

        # Reference table from tests/testthat/test-tidy-layout.R
        expected = {
            1: (6.00, 0.40),
            2: (7.00, 0.50),
            3: (8.00, 0.45),
            4: (9.00, 0.55),
            5: (1.00, 0.45),
            6: (2.00, 0.35),
            7: (0.00, 0.60),
            8: (3.00, 0.50),
            9: (4.00, 0.60),
            10: (5.00, 0.45),
            11: (4.00, 0.00),
            12: (7.25, 0.30),
            13: (8.50, 0.40),
            14: (0.75, 0.20),
            15: (1.50, 0.25),
            16: (4.00, 0.35),
        }
        for i in range(len(tidy_df)):
            r = tidy_df.row_dict(i)
            nid = int(r["node"])
            exp_y, exp_x = expected[nid]
            self.assertAlmostEqual(float(r["y"]), exp_y, places=6)
            self.assertAlmostEqual(float(r["x"]), exp_x, places=6)
            self.assertFalse(math.isnan(float(r["branch"])))
            self.assertFalse(math.isnan(float(r["angle"])))

        # Deep unbalanced tree (2000 tips) does not overflow recursion stack
        deep_tr = stree(2000, type="left")
        deep_df = fortify(deep_tr, layout="tidy", ladderize=False)
        self.assertEqual(len(deep_df), 3999)

        # Zero / negative branch lengths
        tr_zero = read_tree(text="((t1:0,t2:0.5):0,(t3:-0.2,t4:0.1):0.3);")
        df_zero = fortify(tr_zero, layout="tidy")
        self.assertTrue(all(math.isfinite(float(v)) for v in df_zero["x"]))
        self.assertTrue(all(math.isfinite(float(v)) for v in df_zero["y"]))

        # Polytomy support
        tr_poly = read_tree(text="((t1:1,t2:1,t3:1):1,(t4:1,t5:1):1);")
        df_poly = fortify(tr_poly, layout="tidy")
        self.assertEqual(len(df_poly), 8)
        self.assertEqual(len(set(round(float(v), 6) for v in df_poly.filter(lambda r: r["isTip"])["y"])), 5)

        # ggtree(tr, layout='tidy') builds cleanly
        p = ggtree(tr, layout="tidy")
        self.assertTrue(is_ggplot(p))
        built = ggplot_build(p)
        self.assertGreater(len(built.data), 0)

    def test_treelayout_equivalence_and_external_callable(self):
        """Port of tests/testthat/test-treelayout-equivalence.R."""
        tr = read_tree(text="(((t1:1,t2:2):3,t3:4):5,(t4:1,t5:1):2);")
        for lay in ("rectangular", "slanted", "circular", "fan", "equal_angle", "daylight", "ape", "tree_and_leaf", "dendrogram", "inward_circular", "radial", "ellipse", "roundrect", "tidy"):
            df = fortify(tr, layout=lay)
            self.assertEqual(len(df), 9)
            self.assertTrue(all(math.isfinite(float(v)) for v in df["x"]))
            self.assertTrue(all(math.isfinite(float(v)) for v in df["y"]))
            p = ggtree(tr, layout=lay)
            self.assertTrue(is_ggplot(p))

        # External layout function
        def custom_layout(phylo):
            n = Ntip_val = len(phylo.tip_label) + phylo.Nnode
            return DataFrame({
                "node": list(range(1, n + 1)),
                "x": [float(i) for i in range(1, n + 1)],
                "y": [float(i * 2) for i in range(1, n + 1)],
            })

        df_custom = fortify(tr, layout=custom_layout)
        self.assertEqual(list(df_custom["x"]), [float(i) for i in range(1, 10)])
        self.assertEqual(list(df_custom["y"]), [float(i * 2) for i in range(1, 10)])

    def test_tanglegram_and_ggdoubletree(self):
        """Port of tests/testthat/test-tanglegram.R."""
        t1 = read_tree(text="((A:1,B:1):1,C:2);")
        t2 = read_tree(text="((A:1,C:1):1,B:2);")
        tg = as_tanglegram(t1, t2)
        p = ggdoubletree(tg) + geom_tanglelink()
        self.assertTrue(is_ggplot(p))
        built = ggplot_build(p)
        # Tanglelink layer should connect 3 shared tips (A, B, C)
        link_df = built.data[-1]
        self.assertEqual(len(link_df), 3)

    def test_collapse_expand_rotate_flip_scale_view_clade(self):
        """Port of test-collapse.R, test-rotate.R, test-scaleClade.R, test-viewClade.R."""
        tr = read_tree(text="((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
        p = ggtree(tr) + geom_tiplab()
        node_de = MRCA(tr, "D", "E")

        # collapse & expand
        p_col = collapse(p, node=node_de, mode="max", fill="steelblue")
        self.assertGreater(len(p_col.layers), len(p.layers))
        p_exp = expand(p_col, node=node_de)
        self.assertTrue(all(not bool(v) for v in p_exp.data.get(".collapsed", [])))

        # scaleClade
        p_sc = scaleClade(p, node=node_de, scale=0.5)
        self.assertTrue(is_ggplot(p_sc))

        # rotate & flip
        c_nodes = child(tr, node_de)
        p_rot = rotate(p, node=node_de)
        p_flp = flip(p, node1=c_nodes[0], node2=c_nodes[1])
        self.assertNotEqual(list(p.data["y"]), list(p_rot.data["y"]))
        self.assertNotEqual(list(p.data["y"]), list(p_flp.data["y"]))

        # viewClade & zoomClade
        p_vc = viewClade(p, node=node_de)
        self.assertIsNotNone(p_vc._xlim)
        self.assertIsNotNone(p_vc._ylim)
        p_zc = zoomClade(p, node=node_de, xexpand=0.2)
        self.assertIsNotNone(p_zc._xlim)

        # extract tip order (`get_taxa_name`)
        taxa = get_taxa_name(p)
        self.assertEqual(len(taxa), 5)
        self.assertEqual(set(taxa), {"A", "B", "C", "D", "E"})

    def test_gheatmap_and_msaplot(self):
        """Port of test-gheatmap.R and test-msaplot.R."""
        tr = read_tree(text="((A:1,B:1):1,(C:1,D:1):1);")
        p = ggtree(tr) + geom_tiplab()
        mat = DataFrame({
            "trait1": [1.0, 2.0, 3.0, 4.0],
            "trait2": [4.0, 3.0, 2.0, 1.0],
        }, rownames=["A", "B", "C", "D"])
        p_hm = gheatmap(p, mat, offset=0.5, width=0.8, colnames=True)
        self.assertTrue(is_ggplot(p_hm))
        pos = get_heatmap_column_position(p_hm, by="bottom")
        self.assertEqual(list(pos["label"]), ["trait1", "trait2"])

        # msaplot with bundled sample FASTA
        fas = system_file("extdata/sample.fasta", package="ggtree")
        tr_pa = read_tree(system_file("extdata/pa.nwk", package="ggtree"))
        p_msa = msaplot(ggtree(tr_pa), fasta=fas, offset=0.1, width=1.2)
        self.assertTrue(is_ggplot(p_msa))

    def test_cladelab_hilight_strip_taxalink_range(self):
        """Port of test-cladelab-rectangular.R, test-geom_hilight.R, test-geom_strip.R, test-geom_taxalink.R, test-geom_range.R."""
        tr = read_tree(text="((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
        mrca_de = MRCA(tr, "D", "E")
        p = (
            ggtree(tr)
            + geom_tiplab()
            + geom_hilight(node=mrca_de, fill="steelblue", alpha=0.4)
            + geom_cladelab(node=mrca_de, label="Clade DE", offset=0.2)
            + geom_strip("A", "B", label="Strip AB", offset=0.2)
            + geom_taxalink("A", "E", color="red")
            + geom_treescale()
            + geom_rootedge(rootedge=0.5)
        )
        built = ggplot_build(p)
        self.assertEqual(len(built.data), len(p.layers))

        # Beast tree with geom_range
        beast_file = system_file("extdata/MCC_FluA_H3.tree", package="ggtree")
        beast_tr = read_beast(beast_file)
        p_beast = ggtree(beast_tr, mrsd="2013-01-01") + geom_range("height_0.95_HPD", color="red", alpha=0.6) + theme_tree2()
        self.assertTrue(is_ggplot(p_beast))
        built_beast = ggplot_build(p_beast)
        self.assertGreater(len(built_beast.data[0]), 50)

    def test_operators_attach_and_replace_tree(self):
        """Test %<+% (via %, <<, @, +) and %<% tree replacement."""
        tr1 = read_tree(text="((A:1,B:1):1,(C:1,D:1):1);")
        tr2 = read_tree(text="((A:2,C:2):1,(B:2,D:2):1);")
        meta = DataFrame({"taxa": ["A", "B", "C", "D"], "score": [10, 20, 30, 40], "group": ["g1", "g1", "g2", "g2"]})

        p = ggtree(tr1) % meta + geom_tiplab(aes(color=group)) + geom_tippoint(aes(size=score))
        self.assertIn("score", p.data.columns)
        self.assertIn("group", p.data.columns)

        # Replace tree with %<% equivalent (`p % tr2`)
        p2 = p % tr2
        self.assertEqual(len(p2.layers), len(p.layers))
        self.assertEqual(len(p2.data), len(p.data))

    def test_svg_and_png_rendering_and_ggsave(self):
        """Test standalone SVG and PNG rendering and ggsave."""
        tr = rtree(12)
        p = (
            ggtree(tr)
            + geom_tiplab()
            + geom_nodepoint(color="steelblue", size=3)
            + geom_tippoint(color="firebrick", size=3)
            + geom_treescale()
            + ggtitle("Test Phylogenetic Tree")
        )
        svg = p.to_svg()
        self.assertTrue(svg.startswith("<svg"))
        self.assertIn("Test Phylogenetic Tree", svg)

        png = p.to_png(width=400, height=300)
        self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))

        with tempfile.TemporaryDirectory() as tmpdir:
            svg_path = Path(tmpdir) / "tree.svg"
            png_path = Path(tmpdir) / "tree.png"
            ggsave(svg_path, p)
            ggsave(png_path, p)
            self.assertTrue(svg_path.exists() and svg_path.stat().st_size > 200)
            self.assertTrue(png_path.exists() and png_path.stat().st_size > 200)

    def test_aes_subset_complete(self):
        """Port of tests/testthat/test-aes-subset.R."""
        set_seed(1)
        tr = rtree(10)
        group_character = ["group1"] * 3 + ["group2"] * 4 + ["group3"] * 3
        tree_info = DataFrame({
            "label": [f"t{i}" for i in range(1, 11)],
            "group_character": group_character,
            "group_factor": factor(group_character),
            "group_integer": [1] * 3 + [2] * 4 + [3] * 3,
            "in_group2": [g == "group2" for g in group_character],
        })
        tr_plus_info = left_join(tr, tree_info, by="label")

        def n_subsetted(p):
            b = ggplot_build(p)
            layers = [d for d in b.data if "subset" in d.columns]
            self.assertEqual(len(layers), 1)
            return len(layers[0])

        # String literals (#705)
        self.assertEqual(n_subsetted(ggtree(tr_plus_info) + geom_tippoint(aes(subset='group_character == "group2"'), size=3)), 4)
        self.assertEqual(n_subsetted(ggtree(tr_plus_info) + geom_tiplab(aes(subset='group_character == "group2"'))), 4)

        # Factor comparison
        self.assertEqual(n_subsetted(ggtree(tr_plus_info) + geom_tippoint(aes(subset='group_factor == "group2"'), size=3)), 4)

        # Logical column & integer comparison
        self.assertEqual(n_subsetted(ggtree(tr_plus_info) + geom_tippoint(aes(subset="in_group2"), size=3)), 4)
        self.assertEqual(n_subsetted(ggtree(tr_plus_info) + geom_tippoint(aes(subset="group_integer == 2"), size=3)), 4)

        # Local caller variable & nested scope (#705)
        group_of_interest = "group2"
        self.assertEqual(n_subsetted(ggtree(tr_plus_info) + geom_tippoint(aes(subset="group_character == group_of_interest"), size=3)), 4)
        self.assertEqual(n_subsetted(ggtree(tr_plus_info) + geom_tiplab(aes(subset="group_character == group_of_interest"))), 4)

        def nested():
            wanted = "group3"
            return n_subsetted(ggtree(tr_plus_info) + geom_tippoint(aes(subset="group_character == wanted"), size=3))

        self.assertEqual(nested(), 3)

        # TreeDataList + facet_wrap(~.id) (#699)
        set_seed(2020)
        x_tr = rtree(30)
        d = DataFrame({"label": x_tr.tip.label, "var1": np.abs(rnorm(30)), "var2": np.abs(rnorm(30))})
        td_full = full_join(x_tr, d, by="label")
        trs = TreeDataList(TREE1=td_full, TREE2=td_full)
        p_multi = (
            ggtree(trs)
            + facet_wrap("~.id")
            + geom_tippoint(aes(subset=(_id == "TREE1"), colour=SymExpr("var1")))
            + scale_colour_gradient(low="blue", high="red")
            + new_scale_colour()
            + geom_tippoint(aes(colour=SymExpr("var2")), data=td_filter(_id == "TREE2"))
            + scale_colour_viridis_c()
        )
        b_multi = ggplot_build(p_multi)
        point_layers = [d_lay for d_lay in b_multi.data if "shape" in d_lay.columns]
        self.assertEqual(len(point_layers), 2)
        self.assertTrue(all(len(pl) == 30 for pl in point_layers))

    def test_collapse_complete(self):
        """Port of tests/testthat/test-collapse.R."""
        def find_collapse_candidate(plot):
            df = plot.data
            internal_nodes = [int(n) for n in df[~df["isTip"]]["node"].tolist()]
            tip_ys_all = [float(r["y"]) for r in df.to_records() if bool(r["isTip"]) and not is_na(r["y"])]
            for nd in internal_nodes:
                sp = offspring(df, nd)
                if len(sp) == 0:
                    continue
                sp_tip_ys = [float(r["y"]) for r in sp.to_records() if bool(r["isTip"]) and not is_na(r["y"])]
                sp_all_ys = [float(v) for v in sp["y"].tolist() if not is_na(v)]
                if not sp_tip_ys or not sp_all_ys:
                    continue
                tip_span = max(sp_tip_ys) - min(sp_tip_ys)
                max_sp_y = max(sp_all_ys)
                if tip_span > 1 and any(ty > max_sp_y for ty in tip_ys_all):
                    return nd, sp, max_sp_y, tip_span
            self.fail("no suitable internal node found for collapse test")

        set_seed(1)
        tree1 = rtree(12)
        p1 = ggtree(tree1) + geom_tiplab()
        cand_node, sp1, max_sp_y1, tip_span1 = find_collapse_candidate(p1)
        upper_nodes = [int(r["node"]) for r in p1.data.to_records() if bool(r["isTip"]) and float(r["y"]) > max_sp_y1]

        p2 = collapse(p1, cand_node, mode="mixed", height=1)
        tri = p2.layers[-1].data
        tri_ys = [float(v) for v in tri["y"].tolist()]
        self.assertAlmostEqual(max(tri_ys) - min(tri_ys), 1.0, places=8)

        p1_map = {int(r["node"]): float(r["y"]) for r in p1.data.to_records() if not is_na(r["y"])}
        p2_map = {int(r["node"]): float(r["y"]) for r in p2.data.to_records() if not is_na(r["y"])}
        for un in upper_nodes:
            self.assertAlmostEqual(p2_map[un], p1_map[un] - (tip_span1 - 1.0), delta=1e-8)

        # rel() scaling
        p2_rel = collapse(p1, cand_node, mode="mixed", height=rel(0.2))
        tri_rel_ys = [float(v) for v in p2_rel.layers[-1].data["y"].tolist()]
        self.assertAlmostEqual(max(tri_rel_ys) - min(tri_rel_ys), tip_span1 * 0.2, delta=1e-8)

        # expand restores original coordinates
        set_seed(2)
        tree2 = rtree(14)
        p_orig = ggtree(tree2) + geom_tiplab()
        cand2, _, _, _ = find_collapse_candidate(p_orig)
        p_col2 = collapse(p_orig, cand2, mode="mixed", height=1.5)
        p_exp2 = expand(p_col2, cand2)
        for i in range(len(p_orig.data)):
            r_orig = p_orig.data.row_dict(i)
            r_exp = p_exp2.data.row_dict(i)
            self.assertEqual(int(r_exp["node"]), int(r_orig["node"]))
            self.assertEqual(int(r_exp["parent"]), int(r_orig["parent"]))
            self.assertEqual(r_exp["label"], r_orig["label"])
            for col in ("x", "y", "branch", "angle"):
                self.assertAlmostEqual(float(r_exp[col]), float(r_orig[col]), delta=1e-8)

        # height validation error
        p_val = ggtree(rtree(8))
        int_node = int(p_val.data[~p_val.data["isTip"]]["node"].tolist()[0])
        expected_msg = "`height` must be NULL, `ggplot2::rel()` with a single non-negative value, or a single non-negative numeric value."
        with self.assertRaisesRegex(ValueError, "non-negative") as ctx1:
            collapse(p_val, int_node, mode="mixed", height=-1)
        self.assertEqual(str(ctx1.exception), expected_msg)
        with self.assertRaisesRegex(ValueError, "non-negative") as ctx2:
            collapse(p_val, int_node, mode="mixed", height=rel(-0.1))
        self.assertEqual(str(ctx2.exception), expected_msg)

    def test_conversion_as_polytomy(self):
        """Port of tests/testthat/test-conversion.R."""
        f = system_file("extdata/RAxML", "RAxML_bipartitions.H3", package="treeio")
        tree = read_tree(f)
        tree2 = as_polytomy(tree, "node.label", lambda val: float(val) < 70)
        self.assertLess(tree2.Nnode, tree.Nnode)
        for lbl in tree2.node_label or []:
            if lbl is not None and str(lbl).strip() != "":
                self.assertGreater(float(lbl), 70.0)

    def test_geom_balance_suite(self):
        """Port of tests/testthat/test-geom_balance.R."""
        tr = read_tree(text="(t2:0.1280947195,(t4:0.2642134866,(t3:0.9758362926,(t1:0.3494729637,t5:0.189841171,t6:0.189841171):0.7222939907):0.3401968146):0.5143072554);")
        p = ggtree(tr) + geom_balance(9)
        self.assertTrue(is_ggplot(p))
        b = ggplot_build(p)
        self.assertEqual(len(b.data[-1]), 2)

    def test_geom_cladelab_and_cladelabel_parse(self):
        """Port of tests/testthat/test-geom_cladelab.R and test-geom_cladelabel.R."""
        set_seed(1)
        tr = rtree(5)

        def get_cladelab_text(p):
            b = ggplot_build(p)
            for lay, d in zip(p.layers, b.data):
                if lay.geom_type == "cladelab_text" and len(d) > 0:
                    return d
            self.fail("cladelab_text layer not found")

        p_parse = ggtree(tr) + geom_cladelab(node=6, label="italic('Test')", parse=True) + xlim(0, 3)
        t_parse = get_cladelab_text(p_parse)
        self.assertEqual(t_parse["parse"].tolist()[0], True)
        self.assertEqual(t_parse["label"].tolist()[0], "italic('Test')")

        p_lit = ggtree(tr) + geom_cladelab(node=6, label="italic('Test')", parse=False) + xlim(0, 3)
        t_lit = get_cladelab_text(p_lit)
        self.assertEqual(t_lit["parse"].tolist()[0], False)
        self.assertEqual(t_lit["label"].tolist()[0], "italic('Test')")

        p_circ = ggtree(tr, layout="circular") + geom_cladelab(node=6, label="italic('Test')", parse=True)
        t_circ = get_cladelab_text(p_circ)
        self.assertEqual(t_circ["parse"].tolist()[0], True)

        # geom_cladelabel
        tr30 = rtree(30)
        p0 = ggtree(tr30)
        nd = int(p0.data[~p0.data["isTip"]]["node"].tolist()[0])
        p_lbl = p0 + geom_cladelabel(node=nd, label='paste(italic("species name"), "accession number")', parse=True)
        self.assertTrue(is_ggplot(p_lbl))
        ggplot_build(p_lbl)

    def test_geom_hilight_complete(self):
        """Port of tests/testthat/test-geom_hilight.R."""
        set_seed(102)
        tree = rtree(60)
        dat = DataFrame({"id": [62, 88], "type": ["A", "B"]})

        p_grad = ggtree(tree) + geom_hilight(data=dat, mapping=aes(node=id, fill=type_), gradient=True, alpha=0.68)
        self.assertTrue(is_ggplot(p_grad))
        self.assertEqual(len(ggplot_build(p_grad).data[-1]), 2)

        p_rr = ggtree(tree) + geom_hilight(data=dat, mapping=aes(node=id, fill=type_), roundrect=True)
        self.assertTrue(is_ggplot(p_rr))
        self.assertEqual(len(ggplot_build(p_rr).data[-1]), 2)

        p_bot = ggtree(tree) + geom_hilight(data=dat, mapping=aes(node=id, fill=type_), to_bottom=True)
        self.assertEqual(p_bot.layers[0].geom_type, "hilight")

    def test_geom_rootedge_complete(self):
        """Port of tests/testthat/test-geom_rootedge.R."""
        set_seed(42)
        phy = rtree(5, rooted=True)
        td = as_treedata(as_tibble(phy))

        d_phylo = fortify(phy, branch_length="none")
        d_td = fortify(td, branch_length="none")
        self.assertIn("branch.length", d_phylo.columns)
        self.assertIn("branch.length", d_td.columns)

        for re_len in (1, 2):
            p = ggtree(td, branch_length="none") + geom_rootedge(rootedge=re_len)
            d = ggplot_build(p).data[-1]
            self.assertEqual(len(d), 1)
            self.assertAlmostEqual(float(d["xend"][0]), float(d["x"][0]) - re_len, places=8)

        # Default rootedge=None on branch.length='none'
        ggplot_build(ggtree(td, branch_length="none") + geom_rootedge())

        re_p = ggplot_build(ggtree(phy, branch_length="none") + geom_rootedge(rootedge=1)).data[-1]
        re_t = ggplot_build(ggtree(td, branch_length="none") + geom_rootedge(rootedge=1)).data[-1]
        self.assertEqual(re_t["x"].tolist(), re_p["x"].tolist())
        self.assertEqual(re_t["xend"].tolist(), re_p["xend"].tolist())

        # Branch length honored when not 'none'
        d_bl = ggplot_build(ggtree(td) + geom_rootedge(rootedge=1)).data[-1]
        self.assertAlmostEqual(float(d_bl["xend"][0]), float(d_bl["x"][0]) - 1.0, places=8)

    def test_geom_strip_complete(self):
        """Port of tests/testthat/test-geom_strip.R."""
        set_seed(123)
        tree = rtree(10, rooted=True, tip_label=[str(i) for i in range(1, 11)])

        def strip_text_angle(p):
            b = ggplot_build(p)
            for lay, d in zip(p.layers, b.data):
                if lay.geom_type == "strip_text":
                    return float(d["angle"][0])
            self.fail("strip_text layer not found")

        for lay in ("circular", "rectangular"):
            p = ggtree(tree, layout=lay, branch_length="none") + geom_strip(1, 2, label="Test Label", angle="auto")
            ang = strip_text_angle(p)
            self.assertTrue(math.isfinite(ang))

        p_up = ggtree(tree, layout="circular", branch_length="none") + geom_strip(1, 2, label="Test Label", angle="auto", horizontal=True)
        p_rot = ggtree(tree, layout="circular", branch_length="none") + geom_strip(1, 2, label="Test Label", angle="auto", horizontal=False)
        self.assertAlmostEqual(strip_text_angle(p_rot), strip_text_angle(p_up) + 270.0, places=6)

        p_num = ggtree(tree, layout="circular", branch_length="none") + geom_strip(1, 2, label="Test Label", angle=45)
        self.assertAlmostEqual(strip_text_angle(p_num), 45.0, places=6)

    def test_geom_tiplab_complete(self):
        """Port of tests/testthat/test-geom_tiplab.R."""
        def leader_line_layer(b):
            layers = [d for d in b.data if "xend" in d.columns and "label" in d.columns]
            self.assertEqual(len(layers), 1)
            return layers[0]

        tree = rtree(10)
        p = ggtree(tree) + geom_tiplab(align=True, linetype="dotted")
        ll1 = leader_line_layer(ggplot_build(p))
        self.assertEqual(len(ll1), 10)

        d = DataFrame({"label": tree.tip.label, "keep": [True, False] * 5})
        tree2 = left_join(tree, d, by="label")
        p2 = ggtree(tree2) + geom_tiplab(aes(subset="keep"), align=True, linetype="dotted")
        ll2 = leader_line_layer(ggplot_build(p2))
        self.assertEqual(len(ll2), 5)
        global_max_x = max(float(v) for v in p2.data[p2.data["isTip"]]["x"].tolist())
        self.assertTrue(all(abs(float(v) - global_max_x) < 1e-8 for v in ll2["xend"].tolist()))

    def test_geom_tree_complete(self):
        """Port of tests/testthat/test-geom_tree.R."""
        tree = rtree(10)
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            for p_item in (
                ggtree(tree, size=0.5),
                ggtree(tree, size=0.5, ladderize=False),
                ggtree(tree) + geom_tree(size=2),
                ggtree(tree, linewidth=0.5),
                ggtree(tree, size=1, linewidth=4),
                ggtree(tree, linewidth=4, size=1),
            ):
                ggplot_build(p_item)
            for lay in ("rectangular", "dendrogram", "circular", "fan", "slanted", "equal_angle", "ellipse"):
                ggplot_build(ggtree(tree, layout=lay, size=0.5))
            self.assertEqual(len(rec), 0)

        self.assertEqual(set(ggplot_build(ggtree(tree, size=2)).data[0]["linewidth"].tolist()), {2.0})
        self.assertEqual(set(ggplot_build(ggtree(tree, linewidth=2)).data[0]["linewidth"].tolist()), {2.0})
        self.assertEqual(set(ggplot_build(ggtree(tree)).data[0]["linewidth"].tolist()), {0.5})
        self.assertEqual(set(ggplot_build(ggtree(tree, size=1, linewidth=4)).data[0]["linewidth"].tolist()), {4.0})
        self.assertEqual(set(ggplot_build(ggtree(tree, linewidth=4, size=1)).data[0]["linewidth"].tolist()), {4.0})

        b_par = ggplot_build(ggtree(tree, colour="red", alpha=0.3, linetype=2))
        self.assertEqual(set(b_par.data[0]["colour"].tolist()), {"red"})
        self.assertEqual(set(b_par.data[0]["alpha"].tolist()), {0.3})
        self.assertEqual(set(b_par.data[0]["linetype"].tolist()), {2})

        b_pt = ggplot_build(ggtree(tree) + geom_tippoint(size=3))
        tip_pts = [d for d in b_pt.data if "shape" in d.columns]
        self.assertEqual(set(tip_pts[0]["size"].tolist()), {3})

        b_rp = ggplot_build(ggtree(tree) + geom_rootpoint())
        rp = [d for d in b_rp.data if "shape" in d.columns][0]
        self.assertEqual(rp["shape"].tolist()[0], 19)
        self.assertEqual(rp["size"].tolist()[0], 1.5)

    def test_gheatmap_complete(self):
        """Port of tests/testthat/test-gheatmap.R."""
        # Warn when row names do not match
        tr4 = rtree(4)
        mat_miss = DataFrame({"a": [1, 2, 3, 4], "b": [5, 6, 7, 8]}, rownames=[f"missing_{i}" for i in range(1, 5)])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            p_miss = gheatmap(ggtree(tr4), mat_miss)
            self.assertTrue(is_ggplot(p_miss))
            self.assertTrue(any("missing from `data` row names" in str(item.message) for item in w))

        # Custom column labels remapped by source column name
        mat_ok = DataFrame({"A": [1, 2, 3, 4], "B": [5, 6, 7, 8]}, rownames=tr4.tip.label)
        p_cust = gheatmap(ggtree(tr4), mat_ok, custom_column_labels={"Renamed B": "B", "Renamed A": "A"})
        self.assertEqual(p_cust.attrs["data_axis"]["custom_labels"], ["Renamed A", "Renamed B"])

        # Mapping does not leak into plot (#686), tree layer stays 9 rows for 5-tip tree
        tr5 = rtree(5)
        mat5 = DataFrame({"A": runif(5), "B": runif(5)}, rownames=tr5.tip.label)
        p5 = gheatmap(ggtree(tr5), mat5)
        self.assertIsInstance(p5.mapping, Aes)
        self.assertEqual(len(ggplot_build(p5).data[0]), 9)

        # Stacked heatmaps via new_scale_fill on circular tree (#700)
        nwk = system_file("extdata", "sample.nwk", package="treeio")
        tree13 = read_tree(nwk)
        circ = ggtree(tree13, layout="circular")
        df1 = DataFrame(
            {
                "first": ["a", "b", "a", "c", "d", "d", "a", "b", "e", "e", "f", "c", "f"],
                "second": ["z", "z", "z", "z", "y", "y", "y", "y", "x", "x", "x", "a", "a"],
            },
            rownames=tree13.tip.label,
        )
        df2 = DataFrame(
            {
                "A": rnorm(13).tolist(),
                "B": rnorm(13).tolist(),
                "C": rnorm(13).tolist(),
            },
            rownames=tree13.tip.label,
        )
        p_s1 = gheatmap(circ, df1, offset=0.8, width=0.2, colnames_angle=95, colnames_offset_y=0.25) + scale_fill_viridis_d(option="D")
        p_s2 = p_s1 + new_scale_fill()
        p_s3 = gheatmap(p_s2, df2, offset=15, width=0.3, colnames_angle=90, colnames_offset_y=0.25) + scale_fill_viridis_c(option="A")
        b_stack = ggplot_build(p_s3)
        self.assertEqual(len(b_stack.data[0]), 25)
        self.assertEqual(len(b_stack.data[1]), 25)

        # Callable returning {'xy': ...}
        p_xy = ggtree(tr5, layout=lambda obj, **kw: {"xy": np.column_stack([np.arange(1, 10), np.arange(1, 10)])})
        self.assertTrue(is_ggplot(p_xy))
        self.assertEqual(get_layout(p_xy), "external")

    def test_layout_attr_and_tip_y_alignment(self):
        """Port of tests/testthat/test-layout-attr.R and tip_y_alignment checks."""
        set_seed(1)
        tr = rtree(8)
        self.assertEqual(get_layout(ggtree(tr)), "rectangular")
        self.assertEqual(get_layout(ggtree(tr, layout="circular")), "circular")
        self.assertEqual(get_layout(ggtree(tr, layout="daylight")), "daylight")
        self.assertEqual(get_layout(ggtree(tr) + layout_dendrogram()), "dendrogram")
        self.assertEqual(get_layout(ggtree(tr) + layout_fan(angle=90)), "fan")
        self.assertEqual(get_layout(ggtree(tr) + layout_circular()), "circular")
        self.assertIsNone(get_layout(DataFrame({"x": [1]})))

        # tip_y_alignment: ok, duplicated, uneven
        set_seed(7)
        tr50 = rtree(50)
        self.assertEqual(tip_y_alignment(ggtree(tr50)), "ok")
        self.assertEqual(tip_y_alignment(ggtree(tr50, layout="tidy", branch_length="none")), "ok")
        self.assertEqual(tip_y_alignment(ggtree(tr50, layout="tidy")), "duplicated")

        dat50 = DataFrame({"V1": runif(50), "V2": runif(50)}, rownames=tr50.tip.label)
        with warnings.catch_warnings(record=True) as w_dup:
            warnings.simplefilter("always")
            gheatmap(ggtree(tr50, layout="tidy"), dat50, width=0.5, colnames=False)
            self.assertTrue(any("share the same vertical position" in str(item.message) for item in w_dup))

        tr3 = read_tree(text="((a:1,b:1):1,c:2);")
        uneven_df = fortify(tr3)
        uneven_df = uneven_df.copy()
        new_ys = [0.0 if int(n) == 1 else (1.0 if int(n) == 2 else (3.0 if int(n) == 3 else float(yv))) for n, yv in zip(uneven_df["node"].tolist(), uneven_df["y"].tolist())]
        uneven_df["y"] = new_ys
        self.assertEqual(tip_y_alignment(uneven_df), "uneven")

    def test_tanglegram_complete(self):
        """Port of tests/testthat/test-tanglegram.R (merging, layer replay, optimization, both-side)."""
        set_seed(1)
        tr1 = rtree(6)
        tr2 = rtree(6)
        tr2.tip.label = sample(tr1.tip.label, len(tr1.tip.label))
        assoc = DataFrame({"left": tr1.tip.label, "right": tr1.tip.label})

        tg = as_tanglegram(tr1, y=tr2, assoc=assoc)
        df = fortify(tg)
        self.assertEqual(sorted(set(df["side"].tolist())), ["left", "right"])
        self.assertEqual(len(set(df["node"].tolist())), len(df))
        links = df.attrs["tangle_assoc"]
        self.assertEqual(len(links), len(assoc))
        for col in ("x", "y", "xend", "yend"):
            self.assertIn(col, links.columns)

        # Replay layers from GGTree inputs
        set_seed(3)
        r1 = rtree(5)
        r2 = rtree(5)
        r2.tip.label = sample(r1.tip.label, len(r1.tip.label))
        assoc5 = DataFrame({"left": r1.tip.label, "right": r2.tip.label})
        p_paired = ggdoubletree(ggtree(r1) + geom_tiplab(), ggtree(r2) + geom_tippoint(), assoc5)
        self.assertGreaterEqual(len(p_paired.layers), 4)
        ggplot_build(p_paired)

        # Optimization reduces crossings to 0
        tr_l = read_tree(text="((a,b),(c,d));")
        tr_r = read_tree(text="((d,c),(b,a));")
        assoc4 = DataFrame({"left": ["a", "b", "c", "d"], "right": ["a", "b", "c", "d"]})
        plain = fortify(as_tanglegram(tr_l, y=tr_r, assoc=assoc4), optimize=False)
        improved = fortify(as_tanglegram(tr_l, y=tr_r, assoc=assoc4), optimize=True)
        self.assertFalse(plain.attrs["tangle_optimize"]["optimized"])
        self.assertTrue(improved.attrs["tangle_optimize"]["optimized"])
        self.assertGreaterEqual(plain.attrs["tangle_optimize"]["before"], improved.attrs["tangle_optimize"]["after"])
        self.assertEqual(improved.attrs["tangle_optimize"]["after"], 0)

        # Both-side optimization non-regressive
        tr_l5 = read_tree(text="(((a,b),c),(d,e));")
        tr_r5 = read_tree(text="((e,d),(c,(b,a)));")
        assoc5b = DataFrame({"left": ["a", "b", "c", "d", "e"], "right": ["a", "b", "c", "d", "e"]})
        right_only = fortify(as_tanglegram(tr_l5, y=tr_r5, assoc=assoc5b), optimize=True, optimize_side="right")
        both_side = fortify(as_tanglegram(tr_l5, y=tr_r5, assoc=assoc5b), optimize=True, optimize_side="both")
        diag_r = right_only.attrs["tangle_optimize"]
        diag_b = both_side.attrs["tangle_optimize"]
        self.assertEqual(diag_b["optimize_side"], "both")
        self.assertGreaterEqual(diag_r["after"], diag_b["after"])
        self.assertGreaterEqual(diag_b["iterations"], 1)
        ggplot_build(ggdoubletree(tr_l5, tr_r5, assoc5b, optimize=True, optimize_side="both"))

    def test_tidy_layout_invariants(self):
        """Port of remaining checks in tests/testthat/test-tidy-layout.R."""
        set_seed(20260914)
        for _ in range(5):
            tr = rtree(50)
            a = fortify(tr, layout="rectangular")
            b = fortify(tr, layout="tidy")
            self.assertEqual(a[a["isTip"]]["node"].tolist(), b[b["isTip"]]["node"].tolist())

        set_seed(99)
        for i in range(1, 13):
            tr = rtree(40)
            if i % 4 == 0:
                tr = di2multi(tr, tol=0.5)
            d = fortify(tr, layout="tidy")
            nodes = d["node"].tolist()
            parents = d["parent"].tolist()
            xs = [float(v) for v in d["x"].tolist()]
            ys = [float(v) for v in d["y"].tolist()]
            idx_map = {int(n): idx for idx, n in enumerate(nodes)}
            pxs = [xs[idx_map[int(p)]] for p in parents]
            for k in range(len(nodes)):
                for j in range(len(nodes)):
                    if j != k and abs(ys[j] - ys[k]) < 1e-9:
                        self.assertFalse(pxs[j] < xs[k] - 1e-9 and xs[j] > xs[k] + 1e-9)

    def test_svg_color_mapping_col_and_stacked_nodebar(self):
        """Verify SVG rendering for mapped branch/tiplab colors, geom_col facets, and stacked nodebar."""
        tr = read_tree(text="((A:1,B:1):1,(C:1,D:1):1);")
        tr_c = groupClade(tr, MRCA(tr, "A", "B"))
        p_clade = (
            ggtree(tr_c, aes(color="group"))
            + scale_color_manual(values={"0": "#112233", "1": "#D95F02"})
        )
        svg_clade = p_clade.to_svg()
        self.assertIn("#D95F02", svg_clade)
        self.assertIn("#112233", svg_clade)

        meta = DataFrame({"taxa": ["A", "B", "C", "D"], "group": ["g1", "g1", "g2", "g2"], "ab": [2.5, 4.0, 1.5, 3.0]})
        p_tip = (
            ggtree(tr) % meta
            + geom_tiplab(aes(color="group"))
            + scale_color_manual(values={"g1": "#2563EB", "g2": "#DC2626"})
            + geom_facet(panel="Abundance", data=meta, geom=geom_col, mapping=aes(x="ab", fill="group"))
            + scale_fill_manual(values={"g1": "#93C5FD", "g2": "#FCA5A5"})
        )
        svg_tip = p_tip.to_svg()
        self.assertIn("#2563EB", svg_tip)
        self.assertIn("#DC2626", svg_tip)
        self.assertIn("#93C5FD", svg_tip)
        self.assertIn("#FCA5A5", svg_tip)

        bar_df = DataFrame({"node": [5, 6], "x1": [0.7, 0.3], "x2": [0.3, 0.7]})
        bars_stack = nodebar(bar_df, cols=["x1", "x2"], color={"x1": "#1D4ED8", "x2": "#F59E0B"}, position="stack")
        bars_dodge = nodebar(bar_df, cols=["x1", "x2"], color={"x1": "#1D4ED8", "x2": "#F59E0B"}, position="dodge")
        svg_stack = inset(ggtree(tr), bars_stack, width=0.1, height=0.1).to_svg()
        svg_dodge = inset(ggtree(tr), bars_dodge, width=0.1, height=0.1).to_svg()
        self.assertIn("#1D4ED8", svg_stack)
        self.assertIn("#F59E0B", svg_stack)
        self.assertNotEqual(svg_stack, svg_dodge)

    def test_extract_clade_balance_size_and_facet_legend(self):
        """Verify extract_clade preserves node_label, 2-color geom_balance, mapped point size scaling, and facet legend rendering."""
        import re

        tr = read_tree(text="((A:1,B:1)95:1,(C:1,D:1)88:1)100;")
        sub = extract_clade(tr, MRCA(tr, "A", "B"))
        self.assertIsNotNone(sub.node_label)
        self.assertIn("95", sub.node_label)

        split = MRCA(tr, "A", "D")
        p_bal = ggtree(tr) + geom_balance(split, fill=["#3498DB", "#E67E22"], color=["#1F618D", "#B9770E"])
        svg_bal = p_bal.to_svg()
        self.assertIn("#3498DB", svg_bal)
        self.assertIn("#E67E22", svg_bal)

        meta = DataFrame({"taxa": ["A", "B", "C", "D"], "score": [12, 24, 36, 48], "group": ["g1", "g1", "g2", "g2"]})
        p_sz = ggtree(tr) % meta + geom_tippoint(aes(size="score", color="group"))
        svg_sz = p_sz.to_svg()
        radii = [float(m) for m in re.findall(r'<circle [^>]*r="([0-9.]+)"', svg_sz)]
        self.assertTrue(len(radii) >= 4)
        self.assertTrue(all(r <= 15.0 for r in radii))

        p_fac = (
            ggtree(tr)
            + geom_tiplab()
            + geom_facet(panel="Trait", data=meta, geom=geom_point, mapping=aes(x="score", color="group"))
            + scale_color_manual(values={"g1": "#11AA22", "g2": "#AA1122"})
        )
        svg_fac = p_fac.to_svg()
        self.assertIn("#11AA22", svg_fac)
        self.assertIn("g1", svg_fac)
        self.assertIn("g2", svg_fac)


if __name__ == "__main__":
    unittest.main()
