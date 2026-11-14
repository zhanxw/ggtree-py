"""
Automated test suite executing all R man/*.Rd example codes ported to Python.
Verifies that every R example in ggtree works analogously in Python.
"""

import os
import tempfile
import unittest
import warnings
from ggtree import *


class TestAllRExamples(unittest.TestCase):
    """Executes all 38 R man/*.Rd example blocks in Python."""

    def test_01_date_utilities(self):
        # Decimal2Date & Date2decimal.Rd
        d = decimal2Date(2016.35)
        self.assertEqual(d, "2016-05-08")
        dec = Date2decimal("2016-05-05")
        self.assertAlmostEqual(dec, 2016.3415, places=3)

    def test_02_collapse_and_expand(self):
        # collapse.Rd & expand.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p = ggtree(tree) + geom_tiplab()
        p2 = collapse(p, node=21)
        p3 = expand(p2, node=21)
        self.assertTrue(is_ggplot(p2))
        self.assertTrue(is_ggplot(p3))

    def test_03_facet_plot_and_utilities(self):
        # facet_plot.Rd, geom_facet.Rd, facet-utilities.Rd, facet_widths.Rd
        tr = rtree(10)
        dd = data.frame(id=tr.tip.label, value=abs(rnorm(10)))
        p = ggtree(tr)
        p_fp = facet_plot(p, "Trait", data=dd, geom=geom_barh, mapping=aes(x=value)) if "geom_barh" in globals() else facet_plot(p, "Trait", data=dd, geom=geom_col, mapping=aes(x=value))
        p_fp = facet_labeller(p_fp, c(Tree="phylogeny", Trait="bar"))
        p_fp = facet_widths(p_fp, c(Tree=0.5))
        self.assertTrue(is_ggplot(p_fp))

    def test_04_flip_and_rotate(self):
        # flip.Rd & rotate.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p = ggtree(tree) + geom_text2(aes(subset=~isTip, label=node))
        p_flip = flip(p, 17, 21)
        p_rot = rotate(p, 21)
        self.assertTrue(is_ggplot(p_flip))
        self.assertTrue(is_ggplot(p_rot))

    def test_05_geom_aline_and_balance(self):
        # geom_aline.Rd & geom_balance.Rd
        y_tr = rtree(10)
        p1 = ggtree(y_tr) + geom_tiplab(align=TRUE)
        p2 = ggtree(y_tr) + geom_tiplab(align=TRUE, linetype=None) + geom_aline(linetype="dashed")
        p3 = ggtree(y_tr) + geom_balance(16)
        self.assertTrue(is_ggplot(p1))
        self.assertTrue(is_ggplot(p2))
        self.assertTrue(is_ggplot(p3))

    def test_06_geom_cladelab_and_cladelabel(self):
        # geom_cladelab.Rd & geom_cladelabel.Rd
        set_seed(2015 - 12 - 21)
        tree = rtree(30)
        data_df = data.frame(id=c(34, 56), annote=c("another clade", "more clade names"), image=c("a", "b"), group=c("A", "B"), offset=c(0.1, 0.1), offset_text=c(0.1, 0.2))
        p = ggtree(tree) + xlim(NA, 6)
        p2 = (
            p
            + geom_cladelab(node=45, label="test label")
            + geom_cladelab(data=data_df, mapping=aes(node=id, label=annote, image=image, color=group, offset=offset), geom="shadowtext", horizontal=FALSE, angle=90, show_legend=TRUE)
        )
        p3 = ggtree(tree) + xlim(NA, 6) + geom_cladelabel(node=45, label="test label") + geom_cladelabel(node=34, label="another clade")
        self.assertTrue(is_ggplot(p2))
        self.assertTrue(is_ggplot(p3))

    def test_07_geom_hilight(self):
        # geom_hilight.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p1 = ggtree(tree) + geom_hilight(21)
        dat = data.frame(id=c(17, 21), type=c("A", "B"))
        p2 = ggtree(tree) + geom_hilight(data=dat, mapping=aes(node=id, fill=type_))
        p3 = ggtree(tree) + geom_hilight(data=dat, mapping=aes(node=id, fill=type_), to_bottom=TRUE)
        p4 = ggtree(tree, layout="circular") + geom_hilight(data=dat, mapping=aes(node=id, fill=type_), type="rect")
        self.assertTrue(all(is_ggplot(x) for x in (p1, p2, p3, p4)))

    def test_08_geom_label2_point2_segment2_text2(self):
        # geom_label2.Rd, geom_point2.Rd, geom_segment2.Rd, geom_text2.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p = (
            ggtree(tree)
            + geom_point2(aes(subset=~isTip), color="steelblue", size=3)
            + geom_label2(aes(label=label, subset=(isTip == TRUE)))
            + geom_text2(aes(subset=~isTip, label=node), hjust=-0.3)
            + geom_segment2(aes(x=x, xend=x + 1, y=y, yend=y, subset=isTip))
        )
        built = ggplot_build(p)
        self.assertEqual(len(built.data), 6)

    def test_09_geom_nodelab_tiplab_rootedge_rootpoint(self):
        # geom_nodelab.Rd, geom_tiplab.Rd, geom_tiplab2.Rd, geom_rootedge.Rd, geom_rootpoint.Rd
        tr = rtree(10)
        tr.node.label = paste0("node", seq(1, tr.Nnode))
        tr.root.edge = 2
        p = (
            ggtree(tr)
            + geom_tiplab()
            + geom_nodelab()
            + geom_rootedge()
            + geom_rootpoint(color="red", size=3)
            + geom_tippoint()
            + geom_nodepoint()
        )
        p_circ = ggtree(tr, layout="circular") + geom_tiplab2()
        self.assertTrue(is_ggplot(p))
        self.assertTrue(is_ggplot(p_circ))

    def test_10_geom_range_and_beast(self):
        # geom_range.Rd
        file = system.file("extdata/MEGA7", "mtCDNA_timetree.nex", package="treeio")
        beast_file = system.file("extdata/MCC_FluA_H3.tree", package="ggtree")
        tr = read.beast(beast_file)
        p = ggtree(tr) + geom_range("height_0.95_HPD", color="red", size=2, alpha=0.5)
        self.assertTrue(is_ggplot(p))

    def test_11_geom_strip_and_taxalink(self):
        # geom_strip.Rd & geom_taxalink.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p1 = ggtree(tree) + geom_tiplab() + geom_strip("E", "G", color="red", label="clade X", offset=1, offset_text=0.5)
        dat = data.frame(ta1=c("A", "B", "C"), ta2=c("J", "K", "M"), group=c("A", "B", "A"))
        p2 = ggtree(tree) + geom_tiplab() + geom_taxalink(data=dat, mapping=aes(taxa1=ta1, taxa2=ta2, color=group))
        p3 = ggtree(tree) + geom_tiplab() + geom_taxalink(taxa1="A", taxa2="E")
        self.assertTrue(all(is_ggplot(x) for x in (p1, p2, p3)))

    def test_12_geom_tree_treescale_zoom_clade(self):
        # geom_tree.Rd, geom_tree2.Rd, geom_treescale.Rd, geom_zoom_clade.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p1 = ggplot(tree) + geom_tree() + theme_tree()
        p2 = ggtree(tree) + geom_treescale()
        p3 = ggtree(tree) + geom_tiplab() + geom_zoom_clade(node=21, xexpand=0.2)
        self.assertTrue(all(is_ggplot(x) for x in (p1, p2, p3)))

    def test_13_get_clade_position_and_taxa_name(self):
        # get_clade_position.Rd & get_taxa_name.Rd
        tree = rtree(30)
        p = ggtree(tree)
        pos = get_clade_position(p, 35)
        self.assertEqual(len(pos), 1)
        taxa = get_taxa_name(p)
        self.assertEqual(len(taxa), 30)

    def test_14_ggdensitree(self):
        # ggdensitree.Rd
        trees = read.tree(system.file("extdata", "ggdensitree_example.tree", package="ggtree"))
        p1 = ggdensitree(trees, alpha=0.3, colour="steelblue") + geom_tiplab(size=3) + hexpand(0.35)
        trees_f = fortify(trees)
        p2 = ggdensitree(trees_f, layout="slanted", alpha=0.3, colour="steelblue") + geom_tiplab(size=1) + hexpand(0.35)
        self.assertTrue(is_ggplot(p1))
        self.assertTrue(is_ggplot(p2))

    def test_15_ggtree_and_all_layouts(self):
        # ggtree.Rd & layout-utilities.Rd & tree_layout.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p_rect = ggtree(tree)
        p_slant = ggtree(tree, layout="slanted")
        p_circ = ggtree(tree, layout="circular")
        p_fan = ggtree(tree, layout="fan", open_angle=120)
        p_inward = ggtree(tree, layout="inward_circular", xlim=15)
        p_dendro = ggtree(tree, layout="dendrogram")
        p_unrooted = ggtree(tree, layout="unrooted")
        p_eq = ggtree(tree, layout="equal_angle")
        p_day = ggtree(tree, layout="daylight")
        p_tidy = ggtree(tree, layout="tidy")

        # layout utilities
        p_c2 = p_rect + layout_circular()
        p_f2 = p_rect + layout_fan(angle=90)
        p_i2 = p_rect + layout_inward_circular(xlim=15)
        p_d2 = p_rect + layout_dendrogram()
        p_r2 = p_circ + layout_rectangular()
        self.assertTrue(all(is_ggplot(x) for x in (p_rect, p_slant, p_circ, p_fan, p_inward, p_dendro, p_unrooted, p_eq, p_day, p_tidy, p_c2, p_f2, p_i2, p_d2, p_r2)))

    def test_16_gheatmap_and_msaplot(self):
        # gheatmap.Rd, get_heatmap_column_position.Rd, msaplot.Rd
        beast_file = system.file("extdata/MCC_FluA_H3.tree", package="ggtree")
        beast_tree = read.beast(beast_file)
        genotype_file = system.file("extdata/Genotype.txt", package="ggtree")
        genotype = read.table(genotype_file, sep="\t", stringsAsFactors=FALSE)
        p = ggtree(beast_tree, mrsd="2013-01-01") + geom_treescale(x=2008, y=1, offset=2)
        p = p + geom_tiplab(size=2)
        p_hm = gheatmap(p, genotype, offset=5, width=0.5, font_size=3, colnames_angle=-45, hjust=0)
        self.assertTrue(is_ggplot(p_hm))

        fasta = system.file("extdata/FluA_H3_AA.fas", package="ggtree")
        p_msa = msaplot(ggtree(beast_tree), fasta)
        self.assertTrue(is_ggplot(p_msa))

    def test_17_inset_nodebar_nodepie(self):
        # inset.Rd, nodebar.Rd, nodepie.Rd
        set_seed(2015 - 12 - 18)
        tr = rtree(15)
        p = ggtree(tr)
        a = runif(14, 0, 0.33)
        b = runif(14, 0, 0.33)
        c_vals = runif(14, 0, 0.33)
        d = 1 - a - b - c_vals
        dat = data.frame(a=a, b=b, c=c_vals, d=d)
        dat["node"] = [15 + i for i in range(1, 15)]
        bars = nodebar(dat, cols=slice(0, 4), position="dodge", color=c(a="blue", b="red", c="green", d="orange"))
        p_in = inset(p, bars, x="branch", width=0.03, height=0.5)
        pies = nodepie(dat, cols=slice(0, 4), color=c(a="blue", b="red", c="green", d="orange"), alpha=0.6)
        p_in2 = p + geom_inset(pies)
        self.assertTrue(is_ggplot(p_in))
        self.assertTrue(is_ggplot(p_in2))

    def test_18_open_tree_rotate_tree_scaleClade_viewClade_zoomClade(self):
        # open_tree.Rd, rotate_tree.Rd, scaleClade.Rd, viewClade.Rd, zoomClade.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p = ggtree(tree) + geom_tiplab()
        p_open = open_tree(p, 180)
        p_rot = rotate_tree(p_open, 90)
        p_scale = scaleClade(p, node=21, scale=0.3)
        p_view = viewClade(p, MRCA(p, tip=c("I", "L")))
        p_zoom = zoomClade(p, node=21, xexpand=0.2)
        self.assertTrue(all(is_ggplot(x) for x in (p_open, p_rot, p_scale, p_view, p_zoom)))

    def test_19_td_filter_td_mutate_td_unnest_and_label_pad_range_format(self):
        # td-filter.Rd, td-mutate.Rd, td-unnest.Rd, label_pad.Rd, range_format.Rd
        tree = rtree(30)
        p = (
            ggtree(tree)
            + geom_point(data=td_filter(~isTip), color="red", size=3)
            + geom_tiplab(aes(label="new_label"), data=td_mutate(new_label=lambda r: label_pad(r["label"], pad=".")))
        )
        self.assertTrue(is_ggplot(p))
        rf = range_format(["[1.2345, 2.3456]"], digits=2)
        self.assertEqual(rf[0], "1.23 - 2.35")

    def test_20_multiplot_and_patchwork_operators(self):
        # multiplot.Rd
        p1 = ggtree(rtree(10))
        p2 = ggtree(rtree(10), layout="circular")
        mp = multiplot(p1, p2, ncol=2)
        self.assertTrue(is_ggplot(mp))
        self.assertTrue(is_ggplot(p1 | p2))
        self.assertTrue(is_ggplot(p1 / p2))
        self.assertTrue(is_ggplot(p1 + p2 + layout_dendrogram()))

    def test_21_as_polytomy_and_subtrees(self):
        # as.polytomy.Rd & subtrees
        nwk = system.file("extdata/RAxML", "RAxML_bipartitions.H3", package="treeio")
        tree = read.tree(nwk)
        tree2 = as_polytomy(tree, "node.label", lambda val: float(val) < 70)
        for lbl in tree2.node.label:
            if lbl is not None and str(lbl).strip() != "":
                self.assertGreater(float(lbl), 70.0)
        tree3 = as_.polytomy(tree, "node.label", lambda val: float(val) < 70)
        self.assertEqual(len(tree3.tip.label), len(tree.tip.label))
        subs = subtrees(rtree(6))
        self.assertGreaterEqual(len(subs), 1)

    def test_22_tanglegram_and_ggdoubletree_and_striplab(self):
        # as_tanglegram.Rd, ggdoubletree.Rd, geom_tanglelink.Rd, geom_striplab.Rd
        set_seed(123)
        tr1 = rtree(5)
        tr2 = rtree(5)
        tr2.tip.label = [tr1.tip.label[i] for i in (0, 2, 4, 1, 3)]
        assoc = data.frame(left=tr1.tip.label, right=tr2.tip.label)

        tg = as_tanglegram(tr1, y=tr2, assoc=assoc)
        df = fortify(tg)
        self.assertEqual(sorted(set(df["side"].tolist())), ["left", "right"])
        p1 = ggdoubletree(tr1, tr2, assoc) + geom_tanglelink()
        p2 = ggdoubletree(tr1, tr2, assoc, optimize=TRUE, optimize_side="both")
        self.assertTrue(is_ggplot(p1))
        self.assertTrue(is_ggplot(p2))
        self.assertGreaterEqual(len(ggplot_build(p1).data), 4)

        tr = rtree(10)
        dt = data.frame(ta1=c("t5", "t1"), ta2=c("t6", "t3"), group=c("A", "B"))
        p_sl = (
            ggtree(tr)
            + geom_tiplab()
            + geom_striplab(data=dt, mapping=aes(taxa1=ta1, taxa2=ta2, label=group, color=group), align=TRUE, show_legend=FALSE)
        )
        self.assertTrue(is_ggplot(p_sl))
        self.assertGreaterEqual(len(ggplot_build(p_sl).data), 4)

    def test_23_revts_ggexpand_xlim_themes(self):
        # revts.Rd, ggexpand.Rd, xlim_tree.Rd, xlim_expand.Rd, theme_tree*.Rd
        tr = rtree(10)
        p = ggtree(tr) + theme_tree2()
        p2 = revts(p) + scale_x_continuous(labels=abs)
        self.assertTrue(is_ggplot(p2))

        x_tr = rtree(20)
        x_tr.tip.label = paste0("Really_Long_Label_", x_tr.tip.label)
        p_exp = ggtree(x_tr) + geom_tiplab() + ggexpand(1.5, side="h")
        self.assertTrue(is_ggplot(p_exp))

        tr30 = rtree(30)
        d = data.frame(label=tr30.tip.label, value=rnorm(30))
        p_dot = (
            ggtree(tr30)
            + geom_tiplab()
            + geom_facet(panel="Dot", data=d, geom=geom_point, mapping=aes(x=value))
            + xlim_tree(6)
            + xlim_expand(c(-10, 10), "Dot")
        )
        self.assertEqual(p_dot.panel_xlims.get("Tree"), (0.0, 6.0))
        self.assertEqual(p_dot.panel_xlims.get("Dot"), (-10.0, 10.0))

        self.assertTrue(is_ggplot(ggtree(tr) + theme_tree()))
        self.assertTrue(is_ggplot(ggtree(tr) + theme_tree2()))
        self.assertTrue(is_ggplot(ggtree(tr, layout="dendrogram") + theme_dendrogram()))
        self.assertTrue(is_ggplot(ggtree(tr) + theme_inset()))

    def test_24_operators_and_reexports(self):
        # ggtree-format-painter.Rd (%<%, %<+%) & reexports.Rd
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        p = ggtree(tree) + geom_tippoint(color="#b5e521", alpha=1 / 4, size=10)
        p_re = p % rtree(30)
        self.assertTrue(is_ggplot(p_re))
        self.assertEqual(len(p_re.data[p_re.data["isTip"]]), 30)

        dd = data.frame(
            taxa=LETTERS[0:13],
            place=c(rep("GZ", 5), rep("HK", 3), rep("CZ", 4), NA),
            value=[round(float(val), 1) for val in abs(rnorm(13, mean=70, sd=10))],
        )
        p_att = ggtree(tree) % dd + geom_text(aes(color=place, label=label), hjust=-0.5)
        self.assertIn("place", p_att.data.columns)
        self.assertTrue(is_ggplot(p_att))

        # Both direct tidytree function call and symbolic usage of isTip / parent
        self.assertTrue(bool(isTip(tree, 1)))
        self.assertGreater(int(parent(tree, 1)), 13)
        p_sym = ggtree(tree) + geom_point2(aes(subset=isTip), size=2)
        self.assertEqual(len(ggplot_build(p_sym).data[-1]), 13)
        m = MRCA(tree, "A", "B")
        self.assertEqual(len(offspring(tree, m)) >= 2, True)
        self.assertGreaterEqual(len(child(tree, m)), 2)

    def test_25_defunct_interactive_and_ggsave(self):
        # gzoom, scale_color_subtree, set_hilight_legend, interactive, ggsave
        nwk = system.file("extdata", "sample.nwk", package="treeio")
        tree = read.tree(nwk)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            gz = gzoom(tree, c("A", "B", "C"))
            self.assertEqual(len(gz), 2)
            self.assertIsNone(scale_color_subtree())
            p_hl = set_hilight_legend(ggtree(tree), c("red", "blue"), c("g1", "g2"))
            self.assertTrue(is_ggplot(p_hl))
            self.assertGreaterEqual(len(w), 3)

        ggtree_set_interactive()
        try:
            p_int = (
                ggtree(tree)
                + geom_point_interactive(aes(tooltip=label), size=2)
                + geom_tiplab_interactive()
                + geom_nodelab_interactive()
                + geom_cladelab_interactive(node=17, label="clade")
            )
            self.assertTrue(is_ggplot(p_int))
            self.assertGreaterEqual(len(ggplot_build(p_int).data), 5)
        finally:
            ggtree_unset_interactive()

        with tempfile.TemporaryDirectory() as tmpdir:
            svg_path = os.path.join(tmpdir, "out.svg")
            png_path = os.path.join(tmpdir, "out.png")
            p_save = ggtree(tree) + geom_tiplab()
            ggsave(svg_path, plot=p_save)
            ggsave(png_path, plot=p_save)
            self.assertTrue(os.path.exists(svg_path) and os.path.getsize(svg_path) > 100)
            self.assertTrue(os.path.exists(png_path) and os.path.getsize(png_path) > 100)


if __name__ == "__main__":
    unittest.main()
