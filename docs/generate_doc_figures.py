"""Generate all 28 documentation figures (dpi=120) into docs/figures/."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import ggtree as gg

OUT_DIR = Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DPI = 120


def save_fig(name: str, plot, width: float, height: float) -> None:
    out_path = OUT_DIR / name
    gg.ggsave(str(out_path), plot=plot, width=width, height=height, dpi=DPI)
    print(f"Saved {out_path.name} ({width * DPI:.0f}x{height * DPI:.0f}px)")


def _make_10tip_tree():
    return gg.read_tree(
        text="((((t1:0.8,t2:1.1):0.9,(t3:0.7,t4:0.9):1.1):1.2,t5:2.4):0.8,"
        "((t6:1.0,t7:0.8):1.1,(t8:1.3,(t9:0.6,t10:0.9):0.7):0.8):0.9);"
    )


def gen_01_all_layouts_grid() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    layouts = [
        ("rectangular", {}),
        ("roundrect", {}),
        ("slanted", {}),
        ("ellipse", {}),
        ("circular", {}),
        ("fan", {"open_angle": 120}),
        ("inward_circular", {"xlim": 85}),
        ("radial", {}),
        ("dendrogram", {}),
        ("equal_angle", {}),
        ("daylight", {}),
        ("ape", {}),
    ]
    panels = []
    for lay, kw in layouts:
        p = gg.ggtree(tr, layout=lay, color="steelblue", **kw)
        if lay in ("circular", "fan", "inward_circular", "radial"):
            p = p + gg.geom_tiplab2(size=2.0) + gg.hexpand(0.15)
        elif lay == "dendrogram":
            p = p + gg.geom_tiplab(size=2.0) + gg.theme_dendrogram()
        else:
            p = p + gg.geom_tiplab(size=2.0) + gg.hexpand(0.2)
        p = p + gg.ggtitle(lay)
        panels.append(p)
    grid = gg.plot_list(*panels, ncol=4)
    save_fig("01_all_layouts_grid.png", grid, width=14.0, height=10.5)


def gen_02_circular_open_rot() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    base = gg.ggtree(tr, layout="circular", color="#2C3E50") + gg.geom_tiplab2(size=2.8) + gg.hexpand(0.15)
    p1 = base + gg.ggtitle("1. Circular (360°)")
    p2 = gg.open_tree(base, 90) + gg.ggtitle("2. open_tree(p, 90)")
    p3 = gg.rotate_tree(gg.open_tree(base, 90), 45) + gg.ggtitle("3. rotate_tree(..., 45)")
    save_fig("02_circular_open_rot.png", p1 | p2 | p3, width=13.5, height=4.6)


def gen_03_dendrogram() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    clade_ae = gg.MRCA(tr, "A", "E")
    clade_fh = gg.MRCA(tr, "F", "H")
    tr_g = gg.groupClade(tr, [clade_ae, clade_fh])
    p_tree = (
        gg.ggtree(tr_g, gg.aes(color="group"), layout="dendrogram")
        + gg.geom_tiplab(size=3.0)
        + gg.scale_color_manual(values={"0": "#34495E", "1": "#E74C3C", "2": "#2980B9"})
        + gg.theme_dendrogram()
        + gg.ggtitle("Phylogenetic Dendrogram (layout='dendrogram')")
    )

    dist = np.array([
        [0.0, 1.8, 5.2, 5.5, 8.0],
        [1.8, 0.0, 5.0, 5.4, 7.8],
        [5.2, 5.0, 0.0, 2.1, 6.4],
        [5.5, 5.4, 2.1, 0.0, 6.6],
        [8.0, 7.8, 6.4, 6.6, 0.0],
    ])
    hc = gg.hclust(dist, method="average", labels=["Gene_1", "Gene_2", "Gene_3", "Gene_4", "Gene_5"])
    p_hc = (
        gg.ggtree(hc, layout="dendrogram", color="#16A085", size=1.0)
        + gg.geom_tiplab(size=3.0)
        + gg.theme_dendrogram()
        + gg.ggtitle("Hierarchical Clustering (hclust)")
    )
    save_fig("03_dendrogram.png", p_tree | p_hc, width=11.0, height=5.0)


def gen_04_annotated_tree() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    c1 = gg.MRCA(tr, "A", "E")
    c2 = gg.MRCA(tr, "F", "H")
    c3 = gg.MRCA(tr, "I", "L")
    p = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_rootedge(rootedge=8.0)
        + gg.geom_hilight(node=c1, fill="#3498DB", alpha=0.25, extend=4.0)
        + gg.geom_hilight(node=c2, fill="#E67E22", alpha=0.25, extend=4.0)
        + gg.geom_tiplab(size=3.2, offset=1.5)
        + gg.geom_tippoint(color="#2980B9", size=2.5)
        + gg.geom_nodepoint(color="#7F8C8D", size=2.0)
        + gg.geom_cladelab(node=c1, label="italic('Clade I')", parse=True, offset=11.0, offset_text=1.8, barsize=1.5, color="#1F618D", fontsize=3.4)
        + gg.geom_cladelab(node=c2, label="bold('Clade II')", parse=True, offset=11.0, offset_text=1.8, barsize=1.5, color="#B9770E", fontsize=3.4)
        + gg.geom_strip(taxa1="I", taxa2="L", label="Clade III", offset=11.0, offset_text=1.8, barsize=1.5, color="#1E8449", fontsize=3.4)
        + gg.geom_taxalink(taxa1="C", taxa2="G", color="#8E44AD", curvature=-0.45, alpha=0.8)
        + gg.geom_treescale(x=5, y=1.5, width=20)
        + gg.hexpand(0.32)
        + gg.ggtitle("Clade Highlights, Clade Labels, Strips & Taxa Links")
    )
    save_fig("04_annotated_tree.png", p, width=9.0, height=6.0)


def gen_05_beast_gheatmap() -> None:
    beast_tr = gg.read_beast(gg.system_file("extdata/MCC_FluA_H3.tree", package="ggtree"))
    genotype = gg.read_table(gg.system_file("extdata/Genotype.txt", package="ggtree"), sep="\t")
    p = (
        gg.ggtree(beast_tr, mrsd="2013-01-01")
        + gg.geom_range("height_0.95_HPD", color="red", size=1.4, alpha=0.5)
        + gg.geom_tiplab(size=2.0)
        + gg.theme_tree2()
    )
    p_hm = gg.gheatmap(
        p,
        genotype,
        offset=6.5,
        width=0.55,
        font_size=2.8,
        colnames_angle=-45,
        hjust=0.0,
        legend_title="Genotype",
    )
    save_fig("05_beast_gheatmap.png", p_hm, width=10.5, height=9.5)


def gen_06_circular_stacked_heatmap() -> None:
    gg.set_seed(2026)
    tr = gg.rtree(20)
    tips = [str(v) for v in gg.fortify(tr)[gg.fortify(tr)["isTip"]]["label"].tolist()]
    cats = ["Alpha", "Beta", "Gamma", "Delta"]
    df_cat = gg.DataFrame(
        {
            "Host": [cats[i % 4] for i in range(len(tips))],
            "Region": [cats[(i * 3 + 1) % 4] for i in range(len(tips))],
        },
        rownames=tips,
    )
    df_num = gg.DataFrame(
        {
            "ExprA": [round(float(np.sin(i * 0.5) * 2.0 + 3.0), 2) for i in range(len(tips))],
            "ExprB": [round(float(np.cos(i * 0.4) * 2.0 + 3.0), 2) for i in range(len(tips))],
        },
        rownames=tips,
    )
    base = gg.ggtree(tr, layout="circular", open_angle=35) + gg.geom_tiplab2(size=2.2)
    p1 = gg.gheatmap(base, df_cat, offset=0.55, width=0.28, font_size=2.5, legend_title="Category")
    p1 = p1 + gg.scale_fill_viridis_d(option="D") + gg.new_scale_fill()
    p2 = gg.gheatmap(p1, df_num, offset=2.05, width=0.28, font_size=2.5, legend_title="Expression")
    p2 = p2 + gg.scale_fill_viridis_c(option="A") + gg.ggtitle("Circular Tree + Dual gheatmap Scales")
    save_fig("06_circular_stacked_heatmap.png", p2, width=9.0, height=7.5)


def gen_07_msaplot() -> None:
    pa_tree = gg.read_tree(gg.system_file("extdata/pa.nwk", package="ggtree"))
    pa_fas = gg.system_file("extdata/pa.fas", package="ggtree")
    p = gg.msaplot(
        gg.ggtree(pa_tree) + gg.geom_tiplab(size=2.8),
        fasta=pa_fas,
        window=(1, 60),
        offset=0.014,
        width=1.35,
    )
    p = p + gg.ggtitle("Multiple Sequence Alignment Alongside Phylogeny (msaplot)")
    save_fig("07_msaplot.png", p, width=10.5, height=5.8)


def gen_08_facet_plot() -> None:
    gg.set_seed(7)
    tr = gg.rtree(12)
    tip_names = gg.get_taxa_name(gg.ggtree(tr))
    trait_df = gg.DataFrame({
        "id": tip_names,
        "value": [round(float(1.5 + 0.6 * i + np.sin(i) * 1.4), 2) for i in range(len(tip_names))],
        "group": ["High" if i % 2 == 0 else "Low" for i in range(len(tip_names))],
    })
    p_facet = (
        gg.ggtree(tr)
        + gg.geom_tiplab(size=3.0)
        + gg.geom_tippoint(color="#2C3E50", size=2.2)
        + gg.geom_facet(panel="Trait", data=trait_df, geom=gg.geom_point, mapping=gg.aes(x="value", color="group"), size=3.2)
        + gg.scale_color_manual(values={"High": "#E74C3C", "Low": "#2980B9"})
        + gg.hexpand(0.22)
        + gg.theme_bw()
    )
    p_facet = gg.facet_widths(p_facet, widths=[1.1, 1.3])
    p_facet = gg.facet_labeller(p_facet, {"Tree": "Phylogeny", "Trait": "Continuous Trait"})
    save_fig("08_facet_plot.png", p_facet, width=9.0, height=5.2)


def gen_09_inset_pies() -> None:
    tr = _make_10tip_tree()
    p_base = gg.ggtree(tr) + gg.geom_tiplab(size=3.0, offset=0.12) + gg.hexpand(0.18)
    df = gg.fortify(tr)
    int_nodes = [int(r["node"]) for r in df.to_records() if not r["isTip"]]
    anc_rows = []
    for idx, nid in enumerate(int_nodes[:5]):
        a = round(0.15 + 0.14 * ((idx + 1) % 5), 2)
        b = round(0.85 - a - 0.15, 2)
        c = round(1.0 - a - b, 2)
        anc_rows.append({"node": nid, "State_A": a, "State_B": b, "State_C": c})
    anc_df = gg.DataFrame(anc_rows)
    pies = gg.nodepie(anc_df, cols=["State_A", "State_B", "State_C"], color={"State_A": "#E74C3C", "State_B": "#3498DB", "State_C": "#2ECC71"}, alpha=0.88)
    p_pie = gg.inset(p_base, pies, width=0.11, height=0.11) + gg.ggtitle("Ancestral State Proportions via nodepie() + inset()")
    save_fig("09_inset_pies.png", p_pie, width=8.2, height=5.5)


def gen_10_collapsed_and_expanded() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    cnode = gg.MRCA(tr, "A", "E")
    tr_g = gg.groupClade(tr, cnode)
    base = (
        gg.ggtree(tr_g, gg.aes(color="group"))
        + gg.geom_tiplab(size=2.6)
        + gg.scale_color_manual(values={"0": "#2C3E50", "1": "#D35400"})
        + gg.hexpand(0.2)
    )
    p_mix = gg.collapse(base, node=cnode, mode="mixed", fill="#3498DB", alpha=0.45, color="#1F618D") + gg.ggtitle("collapse(..., mode='mixed')")
    p_max = gg.collapse(base, node=cnode, mode="max", fill="#E74C3C", alpha=0.45, color="#922B21") + gg.ggtitle("collapse(..., mode='max')")
    p_min = gg.collapse(base, node=cnode, mode="min", fill="#2ECC71", alpha=0.45, color="#1E8449") + gg.ggtitle("collapse(..., mode='min')")
    p_exp = gg.expand(p_mix, node=cnode) + gg.ggtitle("expand(p_mix, node)")
    save_fig("10_collapsed_and_expanded.png", (p_mix | p_max) / (p_min | p_exp), width=11.0, height=8.2)


def gen_11_ggdensitree() -> None:
    trees = gg.read_tree(gg.system_file("extdata", "ggdensitree_example.tree", package="ggtree"))
    p = (
        gg.ggdensitree(trees, alpha=0.28, colour="steelblue")
        + gg.geom_tiplab(size=2.5, offset=0.05)
        + gg.hexpand(0.32)
        + gg.ggtitle("Multi-Tree Ensemble Overlay (ggdensitree)")
    )
    save_fig("11_ggdensitree.png", p, width=8.5, height=6.5)


def gen_13_multiplot() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    p1 = gg.ggtree(tr, color="#2980B9") + gg.geom_tiplab(size=2.6) + gg.hexpand(0.18) + gg.ggtitle("Rectangular")
    p2 = gg.ggtree(tr, layout="slanted", color="#16A085") + gg.geom_tiplab(size=2.6) + gg.hexpand(0.18) + gg.ggtitle("Slanted")
    p3 = gg.ggtree(tr, layout="circular", color="#8E44AD") + gg.geom_tiplab2(size=2.5) + gg.hexpand(0.15) + gg.ggtitle("Circular")
    mp = gg.multiplot(p1, p2, p3, ncol=3, labels=["A", "B", "C"])
    save_fig("13_multiplot.png", mp, width=13.0, height=4.6)


def gen_14_patchwork() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    c1 = gg.MRCA(tr, "A", "E")
    p_rect = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_hilight(node=c1, fill="gold", alpha=0.35, extend=3.0)
        + gg.geom_tiplab(size=2.6)
        + gg.hexpand(0.2)
        + gg.ggtitle("A. Clade Highlight")
    )
    p_circ = (
        gg.open_tree(gg.ggtree(tr, layout="circular", color="#8E44AD") + gg.geom_tiplab2(size=2.4), 70)
        + gg.hexpand(0.15)
        + gg.ggtitle("B. Fan View")
    )
    tips = gg.get_taxa_name(gg.ggtree(tr))
    hm_df = gg.DataFrame(
        {
            "Trait_1": [round(0.5 + 0.1 * i, 2) for i in range(len(tips))],
            "Trait_2": [round(1.8 - 0.12 * i, 2) for i in range(len(tips))],
            "Trait_3": [round(0.4 + 0.4 * ((i * 2) % 4), 2) for i in range(len(tips))],
        },
        rownames=tips,
    )
    p_bot = gg.gheatmap(
        gg.ggtree(tr) + gg.geom_tiplab(size=2.6) + gg.theme_tree2(),
        hm_df,
        offset=10.0,
        width=0.35,
        font_size=2.8,
        legend_title="Score",
    ) + gg.scale_fill_viridis_c(option="D") + gg.ggtitle("C. Tree + Heatmap Composition: (p_rect | p_circ) / p_bot")
    save_fig("14_patchwork.png", (p_rect | p_circ) / p_bot, width=10.5, height=8.4)


def gen_gs_01_basic_trees() -> None:
    tr = gg.read_tree(text="((A:1.0,B:2.0):1.0,(C:1.5,D:0.8):1.2);")
    p1 = gg.ggtree(tr, color="#2C3E50", size=0.9) + gg.geom_tiplab(size=3.4, offset=0.1) + gg.geom_tippoint(color="#2980B9", size=2.8) + gg.theme_tree2() + gg.hexpand(0.25) + gg.ggtitle("Phylogram (read_tree)")
    p2 = gg.ggtree(tr, branch_length="none", color="#27AE60", size=0.9) + gg.geom_tiplab(size=3.4, offset=0.1) + gg.geom_tippoint(color="#27AE60", size=2.8) + gg.hexpand(0.25) + gg.ggtitle("Cladogram (branch_length='none')")
    dmat = np.array([
        [0.0, 1.5, 5.0, 5.2],
        [1.5, 0.0, 4.8, 5.1],
        [5.0, 4.8, 0.0, 2.0],
        [5.2, 5.1, 2.0, 0.0],
    ])
    hc = gg.hclust(dmat, method="average", labels=["S1", "S2", "S3", "S4"])
    p3 = gg.ggtree(hc, color="#8E44AD", size=0.9) + gg.geom_tiplab(size=3.4, offset=0.1) + gg.geom_tippoint(color="#8E44AD", size=2.8) + gg.theme_tree2() + gg.hexpand(0.28) + gg.ggtitle("Hierarchical (hclust)")
    save_fig("gs_01_basic_trees.png", p1 | p2 | p3, width=13.0, height=4.2)


def gen_gs_02_metadata_operator() -> None:
    tr = gg.read_tree(text="(((A:1.2,B:0.9):1.0,(C:1.5,D:1.1):0.8):0.7,(E:1.3,F:1.6):1.1);")
    meta = gg.DataFrame({
        "taxa": ["A", "B", "C", "D", "E", "F"],
        "score": [12, 24, 18, 36, 28, 45],
        "group": ["Alpha", "Alpha", "Beta", "Beta", "Gamma", "Gamma"],
    })
    p = (
        gg.ggtree(tr, color="#4A5568")
        % meta
        + gg.geom_tiplab(gg.aes(color="group"), size=3.4, offset=0.25)
        + gg.geom_tippoint(gg.aes(size="score", color="group"))
        + gg.scale_color_brewer(palette="Set1")
        + gg.hexpand(0.32)
        + gg.ggtitle("Attaching External Metadata via % Operator")
    )
    save_fig("gs_02_metadata_operator.png", p, width=7.5, height=4.5)


def gen_gs_03_td_filter_mutate() -> None:
    tr = gg.read_tree(text="((Alpha:1.0,Beta_sp:1.4):0.8,(Gamma:1.1,Delta_long:1.8):0.9);")
    p = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_point(data=gg.td_filter("isTip"), color="red", size=3.4)
        + gg.geom_point(data=gg.td_filter("!isTip"), color="steelblue", size=2.8)
        + gg.geom_tiplab(
            data=gg.td_mutate(padded=lambda d: gg.label_pad(d["label"], justify="right", pad=".")),
            mapping=gg.aes(label="padded"),
            size=3.2,
            offset=0.15,
        )
        + gg.hexpand(0.45)
        + gg.ggtitle("td_filter() Selection & td_mutate(label_pad())")
    )
    save_fig("gs_03_td_filter_mutate.png", p, width=7.5, height=4.2)


def gen_gs_04_aes_symbolic_subset() -> None:
    tr = gg.read_tree(text="(((A:1.2,B:0.8):1.1,C:1.8):0.9,(D:1.5,(E:0.7,F:1.1):1.0):0.8);")
    p = (
        gg.ggtree(tr, color="#34495E")
        + gg.geom_nodepoint(gg.aes(subset=~gg.isTip), color="orange", size=3.6)
        + gg.geom_text2(gg.aes(subset=~gg.isTip, label=gg.node), hjust=-0.65, color="#D35400", size=3.2)
        + gg.geom_tippoint(color="#2980B9", size=3.0)
        + gg.geom_tiplab(offset=0.22, size=3.3)
        + gg.geom_text2(gg.aes(x=gg.branch, label="branch.length", subset="branch.length > 0"), vjust=-0.6, color="#7F8C8D", size=2.6)
        + gg.hexpand(0.25)
        + gg.ggtitle("Symbolic Aesthetics (~gg.isTip, label=gg.node, x=gg.branch)")
    )
    save_fig("gs_04_aes_symbolic_subset.png", p, width=8.0, height=4.6)


def gen_gs_05_themes_comparison() -> None:
    tr = gg.read_tree(text="(((t1:1.1,t2:0.9):1.0,t3:1.8):0.8,(t4:1.4,(t5:0.7,t6:1.2):0.9):1.0);")
    base = gg.ggtree(tr, color="#2C3E50") + gg.geom_tippoint(color="#2980B9", size=2.6) + gg.geom_tiplab(size=2.8, offset=0.1) + gg.hexpand(0.22)
    p1 = base + gg.theme_tree() + gg.ggtitle("theme_tree()")
    p2 = base + gg.theme_tree2() + gg.ggtitle("theme_tree2()")
    p3 = base + gg.theme_bw() + gg.ggtitle("theme_bw()")
    p4 = base + gg.theme_minimal() + gg.ggtitle("theme_minimal()")
    save_fig("gs_05_themes_comparison.png", (p1 | p2) / (p3 | p4), width=10.0, height=7.0)


def gen_lm_04_mrsd_and_revts() -> None:
    beast_tr = gg.read_beast(gg.system_file("extdata/MCC_FluA_H3.tree", package="ggtree"))
    p_cal = (
        gg.ggtree(beast_tr, mrsd="2013-01-01", color="#2C3E50")
        + gg.geom_tiplab(size=1.8)
        + gg.theme_tree2()
        + gg.hexpand(0.28)
        + gg.ggtitle("Calendar Axis: mrsd='2013-01-01'")
    )
    p_rev = (
        gg.revts(gg.ggtree(beast_tr, color="#1F618D"))
        + gg.geom_tiplab(size=1.8)
        + gg.theme_tree2()
        + gg.scale_x_continuous(labels=abs)
        + gg.hexpand(0.28)
        + gg.ggtitle("Reverse Timescale: revts() + abs Labels")
    )
    save_fig("lm_04_mrsd_and_revts.png", p_cal | p_rev, width=13.5, height=7.5)


def gen_lm_05_group_clade_otu() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    c1 = gg.MRCA(tr, "A", "E")
    c2 = gg.MRCA(tr, "I", "L")
    tr_cl = gg.groupClade(tr, [c1, c2])
    p_cl = (
        gg.ggtree(tr_cl, gg.aes(color="group"), size=1.0)
        + gg.geom_tiplab(gg.aes(color="group"), size=3.0)
        + gg.scale_color_brewer(palette="Set1")
        + gg.hexpand(0.22)
        + gg.ggtitle("groupClade(tr, [c1, c2])")
    )
    tr_otu = gg.groupOTU(tr, {"Clade_X": ["A", "B", "D", "E"], "Clade_Y": ["G", "H", "K"]})
    p_otu = (
        gg.ggtree(tr_otu, gg.aes(color="group"), size=1.0)
        + gg.geom_tiplab(gg.aes(color="group"), size=3.0)
        + gg.scale_color_brewer(palette="Dark2")
        + gg.hexpand(0.22)
        + gg.ggtitle("groupOTU(tr, dict_of_taxa)")
    )
    save_fig("lm_05_group_clade_otu.png", p_cl | p_otu, width=11.5, height=5.2)


def gen_lm_07_rotate_flip_scale() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    c_ae = gg.MRCA(tr, "A", "E")
    c_fh = gg.MRCA(tr, "F", "H")
    base = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_hilight(node=c_ae, fill="#3498DB", alpha=0.25, extend=3.0)
        + gg.geom_hilight(node=c_fh, fill="#E67E22", alpha=0.25, extend=3.0)
        + gg.geom_tiplab(size=2.6)
        + gg.hexpand(0.18)
    )
    p1 = base + gg.ggtitle("1. Original")
    p2 = gg.rotate(base, c_ae) + gg.ggtitle("2. rotate(p, node_AE)")
    p3 = gg.flip(base, c_ae, c_fh) + gg.ggtitle("3. flip(p, AE, FH)")
    p4 = gg.scaleClade(base, c_ae, scale=0.35) + gg.ggtitle("4. scaleClade(p, AE, 0.35)")
    save_fig("lm_07_rotate_flip_scale.png", (p1 | p2) / (p3 | p4), width=11.0, height=8.2)


def gen_lm_08_viewclade_gzoom() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    focus = gg.MRCA(tr, "A", "E")
    p_full = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_hilight(node=focus, fill="steelblue", alpha=0.3, extend=2.0)
        + gg.geom_tiplab(size=2.6)
        + gg.hexpand(0.2)
        + gg.ggtitle("Overview (Highlighted Clade)")
    )
    p_view = (
        gg.viewClade(
            gg.ggtree(tr, color="#1F618D") + gg.geom_tiplab(size=3.0) + gg.hexpand(0.15),
            node=focus,
            xmax_adjust=8.0,
        )
        + gg.ggtitle("viewClade(p, node)")
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        gz = gg.gzoom(tr, focus=focus, subtree=False)
    p_gz_left = gz.plots[0] + gg.ggtitle("gzoom() Overview")
    p_gz_right = gz.plots[1] + gg.ggtitle("gzoom() Focal Clade")
    save_fig("lm_08_viewclade_gzoom.png", (p_full | p_view) / (p_gz_left | p_gz_right), width=11.0, height=8.0)


def gen_lm_09_as_polytomy() -> None:
    raxml_full = gg.read_tree(gg.system_file("extdata", "RAxML", "RAxML_bipartitions.H3", package="treeio"))
    sub_tr = gg.extract_clade(raxml_full, 94)
    poly_tr = gg.as_polytomy(sub_tr, feature="node.label", fun=lambda x: float(x) < 70)
    p_orig = (
        gg.ggtree(sub_tr, color="#2C3E50")
        + gg.geom_nodelab(size=2.6, hjust=-0.15, color="#C0392B")
        + gg.geom_tiplab(size=2.3, offset=0.0005)
        + gg.hexpand(0.45)
        + gg.ggtitle("Original Bifurcating RAxML Subtree")
    )
    p_poly = (
        gg.ggtree(poly_tr, color="#1F618D")
        + gg.geom_nodelab(size=2.6, hjust=-0.15, color="#1E8449")
        + gg.geom_tiplab(size=2.3, offset=0.0005)
        + gg.hexpand(0.45)
        + gg.ggtitle("as_polytomy(tr, fun=lambda x: float(x) < 70)")
    )
    save_fig("lm_09_as_polytomy.png", p_orig | p_poly, width=13.5, height=5.8)


def gen_ac_02_hilight_and_balance() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    c1 = gg.MRCA(tr, "A", "E")
    c2 = gg.MRCA(tr, "F", "H")
    p_hilite = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_hilight(node=c1, fill="#3498DB", alpha=0.35, extend=4.0)
        + gg.geom_hilight(node=c2, fill="#E74C3C", alpha=0.35, extend=4.0)
        + gg.geom_cladelab(node=c1, label="Clade A-E", offset=9.0, barsize=1.5, color="#1F618D", fontsize=3.2)
        + gg.geom_cladelab(node=c2, label="Clade F-H", offset=9.0, barsize=1.5, color="#922B21", fontsize=3.2)
        + gg.geom_tiplab(size=2.8, offset=1.2)
        + gg.hexpand(0.28)
        + gg.ggtitle("geom_hilight() + geom_cladelab()")
    )
    split_node = gg.MRCA(tr, "A", "H")
    p_bal = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_balance(node=split_node, fill=["#3498DB", "#E67E22"], color=["#1F618D", "#B9770E"], alpha=0.32, extend=4.0)
        + gg.geom_nodepoint(gg.aes(subset=f"node == {split_node}"), color="#D35400", size=4.0)
        + gg.geom_tiplab(size=2.8, offset=1.2)
        + gg.hexpand(0.22)
        + gg.ggtitle("geom_balance(node=MRCA(A, H))")
    )
    save_fig("ac_02_hilight_and_balance.png", p_hilite | p_bal, width=11.5, height=5.5)


def gen_ac_03_taxalink_and_strip() -> None:
    tr = gg.read_tree(gg.system_file("extdata", "sample.nwk", package="treeio"))
    p_rect = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_tiplab(size=2.8, offset=1.2)
        + gg.geom_strip(taxa1="A", taxa2="E", label="Strip A-E", offset=10.0, offset_text=1.5, barsize=1.6, color="#2980B9", fontsize=3.2)
        + gg.geom_strip(taxa1="I", taxa2="L", label="Strip I-L", offset=10.0, offset_text=1.5, barsize=1.6, color="#27AE60", fontsize=3.2)
        + gg.geom_taxalink(taxa1="B", taxa2="G", color="#E74C3C", curvature=-0.4, alpha=0.85)
        + gg.geom_taxalink(taxa1="D", taxa2="K", color="#8E44AD", curvature=-0.4, alpha=0.85)
        + gg.hexpand(0.30)
        + gg.ggtitle("Rectangular geom_strip() & geom_taxalink()")
    )
    p_circ = (
        gg.ggtree(tr, layout="circular", open_angle=45, color="#2C3E50")
        + gg.geom_tiplab2(size=2.6)
        + gg.geom_taxalink(taxa1="B", taxa2="G", color="#E74C3C", curvature=0.5, alpha=0.85)
        + gg.geom_taxalink(taxa1="D", taxa2="K", color="#8E44AD", curvature=0.5, alpha=0.85)
        + gg.hexpand(0.18)
        + gg.ggtitle("Circular Taxa Links")
    )
    save_fig("ac_03_taxalink_and_strip.png", p_rect | p_circ, width=11.5, height=5.4)


def gen_ac_07_multifacet_panels() -> None:
    tr = _make_10tip_tree()
    tips = gg.get_taxa_name(gg.ggtree(tr))
    ab_df = gg.DataFrame({
        "id": tips,
        "abundance": [12, 28, 19, 34, 25, 42, 31, 16, 38, 22],
        "habitat": ["Marine", "Marine", "Fresh", "Fresh", "Soil", "Marine", "Soil", "Fresh", "Soil", "Marine"],
    })
    tr_df = gg.DataFrame({
        "id": tips,
        "trait_x0": [0.0] * len(tips),
        "trait_val": [round(1.2 + 0.45 * i + np.cos(i) * 0.8, 2) for i in range(len(tips))],
        "habitat": ["Marine", "Marine", "Fresh", "Fresh", "Soil", "Marine", "Soil", "Fresh", "Soil", "Marine"],
    })
    p = (
        gg.ggtree(tr, color="#2C3E50")
        + gg.geom_tiplab(size=2.9)
        + gg.geom_facet(panel="Abundance", data=ab_df, geom=gg.geom_col, mapping=gg.aes(x="abundance", fill="habitat"))
        + gg.geom_facet(panel="Trait", data=tr_df, geom=gg.geom_segment, mapping=gg.aes(x="trait_x0", xend="trait_val", y="y", yend="y"), color="#7F8C8D")
        + gg.geom_facet(panel="Trait", data=tr_df, geom=gg.geom_point, mapping=gg.aes(x="trait_val", color="habitat"), size=3.2)
        + gg.scale_color_manual(values={"Marine": "#3498DB", "Fresh": "#2ECC71", "Soil": "#E67E22"})
        + gg.scale_fill_manual(values={"Marine": "#3498DB", "Fresh": "#2ECC71", "Soil": "#E67E22"})
        + gg.theme_bw()
        + gg.hexpand(0.2)
    )
    p = gg.facet_widths(p, widths=[1.1, 0.9, 0.9])
    p = gg.facet_labeller(p, {"Tree": "Phylogeny", "Abundance": "Abundance (geom_col)", "Trait": "Trait (segment + point)"})
    save_fig("ac_07_multifacet_panels.png", p, width=11.5, height=5.4)


def gen_ac_09_inset_bars() -> None:
    tr = _make_10tip_tree()
    p_base = gg.ggtree(tr, color="#2C3E50") + gg.geom_tiplab(size=2.8, offset=0.12) + gg.hexpand(0.2)
    df = gg.fortify(tr)
    int_nodes = [int(r["node"]) for r in df.to_records() if not r["isTip"]]
    bar_rows = []
    for idx, nid in enumerate(int_nodes[:5]):
        a = round(0.2 + 0.12 * ((idx + 1) % 4), 2)
        b = round(0.35 + 0.1 * (idx % 3), 2)
        c = round(max(0.1, 1.0 - a - b), 2)
        bar_rows.append({"node": nid, "Env_A": a, "Env_B": b, "Env_C": c})
    bar_df = gg.DataFrame(bar_rows)
    colors = {"Env_A": "#E74C3C", "Env_B": "#3498DB", "Env_C": "#2ECC71"}
    bars_stack = gg.nodebar(bar_df, cols=["Env_A", "Env_B", "Env_C"], color=colors, alpha=0.88, position="stack")
    bars_dodge = gg.nodebar(bar_df, cols=["Env_A", "Env_B", "Env_C"], color=colors, alpha=0.88, position="dodge")
    p_stack = gg.inset(p_base, bars_stack, width=0.09, height=0.10) + gg.ggtitle("nodebar(..., position='stack')")
    p_dodge = gg.inset(p_base, bars_dodge, width=0.11, height=0.10) + gg.ggtitle("nodebar(..., position='dodge')")
    save_fig("ac_09_inset_bars.png", p_stack | p_dodge, width=11.5, height=5.2)


def gen_ac_11_tanglegram_before_after() -> None:
    t1 = gg.read_tree(text="(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);")
    t2 = gg.read_tree(text="(((H:1,C:1):1,(A:1,F:1):1):1,((D:1,G:1):1,(B:1,E:1):1):1);")
    tips = ["A", "B", "C", "D", "E", "F", "G", "H"]
    assoc = gg.DataFrame({"left": tips, "right": tips})
    p1 = gg.ggtree(t1, color="#2980B9") + gg.geom_tiplab(size=3.0)
    p2 = gg.ggtree(t2, color="#16A085") + gg.geom_tiplab(size=3.0)
    p_raw = gg.ggdoubletree(p1, p2, assoc=assoc, gap=0.22, optimize=False, link_color="#E74C3C") + gg.ggtitle("Before: optimize=False (Crossed Links)")
    p_opt = gg.ggdoubletree(p1, p2, assoc=assoc, gap=0.22, optimize=True, optimize_side="both", link_color="#2980B9") + gg.ggtitle("After: optimize=True, optimize_side='both'")
    save_fig("ac_11_tanglegram_before_after.png", p_raw | p_opt, width=12.5, height=5.2)


def main() -> None:
    generators = [
        gen_01_all_layouts_grid,
        gen_02_circular_open_rot,
        gen_03_dendrogram,
        gen_04_annotated_tree,
        gen_05_beast_gheatmap,
        gen_06_circular_stacked_heatmap,
        gen_07_msaplot,
        gen_08_facet_plot,
        gen_09_inset_pies,
        gen_10_collapsed_and_expanded,
        gen_11_ggdensitree,
        gen_13_multiplot,
        gen_14_patchwork,
        gen_gs_01_basic_trees,
        gen_gs_02_metadata_operator,
        gen_gs_03_td_filter_mutate,
        gen_gs_04_aes_symbolic_subset,
        gen_gs_05_themes_comparison,
        gen_lm_04_mrsd_and_revts,
        gen_lm_05_group_clade_otu,
        gen_lm_07_rotate_flip_scale,
        gen_lm_08_viewclade_gzoom,
        gen_lm_09_as_polytomy,
        gen_ac_02_hilight_and_balance,
        gen_ac_03_taxalink_and_strip,
        gen_ac_07_multifacet_panels,
        gen_ac_09_inset_bars,
        gen_ac_11_tanglegram_before_after,
    ]
    for fn in generators:
        fn()
    print(f"Done: {len(generators)} figures saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
