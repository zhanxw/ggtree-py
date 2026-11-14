"""
End-to-end runnable showcase demonstrating R ggtree workflows ported to Python.

Run from the repository root with:
    python3 -s examples/run_all_examples.py
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

# Ensure repository root is importable when invoked as `python3 examples/run_all_examples.py`
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np

from ggtree import (
    DataFrame,
    MRCA,
    NA,
    TRUE,
    aes,
    as_polytomy,
    as_tanglegram,
    c,
    child,
    collapse,
    expand,
    facet_labeller,
    facet_plot,
    facet_widths,
    flip,
    fortify,
    geom_aline,
    geom_balance,
    geom_cladelab,
    geom_col,
    geom_facet,
    geom_hilight,
    geom_inset,
    geom_label2,
    geom_nodelab,
    geom_nodepoint,
    geom_point,
    geom_point2,
    geom_range,
    geom_rootedge,
    geom_rootpoint,
    geom_segment2,
    geom_strip,
    geom_tanglelink,
    geom_taxalink,
    geom_text2,
    geom_tiplab,
    geom_tiplab2,
    geom_tippoint,
    geom_tree,
    geom_treescale,
    get_layout,
    ggdensitree,
    ggdoubletree,
    ggplot,
    ggplot_build,
    ggsave,
    ggtitle,
    ggtree,
    gheatmap,
    groupClade,
    groupOTU,
    hexpand,
    inset,
    is_ggplot,
    label,
    label_pad,
    layout_circular,
    layout_dendrogram,
    layout_fan,
    msaplot,
    multiplot,
    new_scale_fill,
    node,
    nodebar,
    nodepie,
    open_tree,
    plot_list,
    read,
    read_beast,
    read_tree,
    rel,
    revts,
    rnorm,
    rotate,
    rotate_tree,
    rtree,
    runif,
    scaleClade,
    scale_fill_viridis_c,
    scale_fill_viridis_d,
    scale_x_continuous,
    set_seed,
    system_file,
    td_filter,
    td_mutate,
    theme_dendrogram,
    theme_tree,
    theme_tree2,
    tip_y_alignment,
    viewClade,
    xlim,
    xlim_expand,
    xlim_tree,
    zoomClade,
    isTip,
)


def _verify_and_save_outputs(
    plots: dict[str, tuple[object, tuple[float, float]]],
    out_dir: Path,
) -> list[tuple[str, Path, Path]]:
    """Save each named plot to SVG and PNG in out_dir and verify valid file headers."""
    out_dir.mkdir(parents=True, exist_ok=True)
    saved: list[tuple[str, Path, Path]] = []
    for name, (p, (w_in, h_in)) in plots.items():
        svg_path = out_dir / f"{name}.svg"
        png_path = out_dir / f"{name}.png"
        ggsave(svg_path, plot=p, width=w_in, height=h_in, dpi=100)
        ggsave(png_path, plot=p, width=w_in, height=h_in, dpi=100)
        svg_bytes = svg_path.read_text(encoding="utf-8")
        png_bytes = png_path.read_bytes()
        assert svg_bytes.lstrip().startswith("<svg"), f"Invalid SVG header for {name}"
        assert len(svg_bytes) > 200, f"Empty SVG output for {name}"
        assert png_bytes.startswith(b"\x89PNG\r\n\x1a\n"), f"Invalid PNG header for {name}"
        assert len(png_bytes) > 200, f"Empty PNG output for {name}"
        saved.append((name, svg_path, png_path))
    return saved


def main() -> int:
    set_seed(2026)
    outputs_to_render: dict[str, tuple[object, tuple[float, float]]] = {}

    # 1. All major tree layouts and layout metadata
    nwk_file = system_file("extdata", "sample.nwk", package="treeio")
    sample_tree = read_tree(nwk_file)
    lay_plots = []
    for lay in (
        "rectangular",
        "slanted",
        "circular",
        "fan",
        "inward_circular",
        "dendrogram",
        "equal_angle",
        "daylight",
        "tidy",
        "ellipse",
        "roundrect",
        "radial",
    ):
        p_lay = ggtree(sample_tree, layout=lay) + ggtitle(lay)
        assert is_ggplot(p_lay)
        assert get_layout(p_lay) == lay
        assert len(ggplot_build(p_lay).data) >= 1
        lay_plots.append(p_lay)

    assert tip_y_alignment(ggtree(sample_tree)) == "ok"
    outputs_to_render["01_all_layouts_grid"] = (plot_list(*lay_plots, ncol=4), (12.0, 8.5))
    outputs_to_render["02_circular_open_rot"] = (
        rotate_tree(open_tree(ggtree(sample_tree, layout="circular") + geom_tiplab2(), 90), 45)
        + ggtitle("Circular Fan (open=90, rotate=45)"),
        (8.0, 7.5),
    )
    outputs_to_render["03_dendrogram"] = (
        ggtree(sample_tree) + layout_dendrogram() + geom_tiplab() + theme_dendrogram() + ggtitle("Dendrogram Layout"),
        (8.5, 6.0),
    )

    # 2. Core annotations, symbolic subsetting, clades, strips, and links
    tr_bal = read_tree(text="(t2:0.13,(t4:0.26,(t3:0.98,(t1:0.35,t5:0.19):0.72):0.34):0.51);")
    mrca_de = MRCA(sample_tree, "I", "L")
    sample_tree.root_edge = 0.4
    p_ann = (
        ggtree(sample_tree)
        + geom_rootedge()
        + geom_rootpoint(color="darkred", size=2.5)
        + geom_tippoint(color="steelblue", size=2)
        + geom_nodepoint(aes(subset=~isTip), color="orange", size=2)
        + geom_text2(aes(subset=~isTip, label=node), hjust=-0.3, size=2.5)
        + geom_tiplab(align=TRUE, linetype="dotted", offset=0.2)
        + geom_hilight(node=mrca_de, fill="steelblue", alpha=0.35, to_bottom=TRUE)
        + geom_cladelab(node=mrca_de, label="italic('Clade IL')", parse=TRUE, offset=2.0)
        + geom_strip("A", "C", label="Strip AC", offset=2.0, angle="auto")
        + geom_taxalink("A", "E", color="firebrick")
        + geom_treescale()
        + theme_tree2()
        + ggtitle("Annotated Sample Tree")
    )
    outputs_to_render["04_annotated_tree"] = (p_ann, (9.5, 6.5))
    assert is_ggplot(ggtree(tr_bal) + geom_balance(8))

    # 3. BEAST time-scaled tree + geom_range + revts + gheatmap + new_scale_fill
    beast_file = system_file("extdata/MCC_FluA_H3.tree", package="ggtree")
    beast_tree = read_beast(beast_file)
    genotype_file = system_file("extdata/Genotype.txt", package="ggtree")
    genotype = read.table(genotype_file, sep="\t", stringsAsFactors=False)

    p_beast = (
        ggtree(beast_tree, mrsd="2013-01-01")
        + geom_range("height_0.95_HPD", color="red", size=1.5, alpha=0.5)
        + geom_tiplab(size=2)
        + theme_tree2()
    )
    p_hm = (
        gheatmap(p_beast, genotype, offset=5, width=0.5, font_size=3, colnames_angle=-45, hjust=0)
        + ggtitle("FluA H3 BEAST Tree + Genotype Heatmap")
    )
    assert is_ggplot(revts(ggtree(beast_tree)) + scale_x_continuous(labels=abs))
    outputs_to_render["05_beast_gheatmap"] = (p_hm, (11.0, 11.0))

    # Stacked circular gheatmaps via new_scale_fill
    df1 = DataFrame(
        {
            "first": ["a", "b", "a", "c", "d", "d", "a", "b", "e", "e", "f", "c", "f"],
            "second": ["z", "z", "z", "z", "y", "y", "y", "y", "x", "x", "x", "a", "a"],
        },
        rownames=sample_tree.tip_label,
    )
    df2 = DataFrame(
        {"M1": rnorm(13).tolist(), "M2": rnorm(13).tolist()},
        rownames=sample_tree.tip_label,
    )
    p_circ_hm = (
        gheatmap(ggtree(sample_tree, layout="circular"), df1, offset=0.8, width=0.2)
        + scale_fill_viridis_d(option="D")
        + new_scale_fill()
    )
    p_circ_hm2 = (
        gheatmap(p_circ_hm, df2, offset=6.0, width=0.3)
        + scale_fill_viridis_c(option="A")
        + ggtitle("Circular Stacked Heatmaps (Viridis + Magma)")
    )
    outputs_to_render["06_circular_stacked_heatmap"] = (p_circ_hm2, (9.0, 8.5))

    # 4. msaplot, facet_plot, and inset charts
    fas_file = system_file("extdata/pa.fas", package="ggtree")
    pa_tree = read_tree(system_file("extdata/pa.nwk", package="ggtree"))
    outputs_to_render["07_msaplot"] = (
        msaplot(ggtree(pa_tree) + geom_tiplab(size=2.5), fasta=fas_file, window=(1, 60), offset=0.1, width=1.2)
        + ggtitle("Multiple Sequence Alignment (PA cols 1..60)"),
        (10.0, 6.5),
    )

    tr10 = rtree(10)
    trait_df = DataFrame({"id": tr10.tip_label, "value": np.abs(rnorm(10)).tolist()})
    p_facet = (
        ggtree(tr10)
        + geom_tiplab()
        + geom_facet(panel="Trait", data=trait_df, geom=geom_point, mapping=aes(x="value"))
        + xlim_tree(6)
        + xlim_expand(c(-3, 3), "Trait")
        + ggtitle("Facet Plot (Phylogeny + Trait)")
    )
    p_facet = facet_widths(facet_labeller(p_facet, c(Tree="Phylogeny", Trait="Trait")), c(Tree=0.6))
    outputs_to_render["08_facet_plot"] = (p_facet, (9.5, 6.0))

    a_vals = runif(9, 0, 0.33)
    b_vals = runif(9, 0, 0.33)
    c_vals = runif(9, 0, 0.33)
    d_vals = 1.0 - a_vals - b_vals - c_vals
    pie_df = DataFrame({"a": a_vals, "b": b_vals, "c": c_vals, "d": d_vals, "node": list(range(11, 20))})
    pies = nodepie(pie_df, cols=range(1, 5), color=c(a="blue", b="red", c="green", d="orange"), alpha=0.7)
    outputs_to_render["09_inset_pies"] = (
        inset(ggtree(tr10) + geom_tiplab(), pies, width=0.06, height=0.06) + ggtitle("Node Pie Chart Insets"),
        (8.5, 6.5),
    )

    # 5. Manipulations, grouping, polytomy collapsing, and operators (% / format painter)
    p_col = collapse(ggtree(sample_tree) + geom_tiplab(), node=mrca_de, mode="mixed", height=rel(0.35), fill="steelblue") + ggtitle("Collapsed (mixed)")
    p_exp = expand(p_col, node=mrca_de) + ggtitle("Expanded")
    p_scaled = scaleClade(ggtree(sample_tree), node=mrca_de, scale=0.4)
    ch = child(sample_tree, mrca_de)
    p_flip = flip(rotate(ggtree(sample_tree), mrca_de), ch[0], ch[1])
    p_view = zoomClade(viewClade(ggtree(sample_tree), mrca_de), mrca_de, xexpand=0.2)
    outputs_to_render["10_collapsed_and_expanded"] = (p_col | p_exp, (10.5, 6.0))
    assert all(is_ggplot(x) for x in (p_scaled, p_flip, p_view))

    raxml_tree = read_tree(system_file("extdata/RAxML", "RAxML_bipartitions.H3", package="treeio"))
    poly_tree = as_polytomy(raxml_tree, "node.label", lambda val: float(val) < 70)
    assert poly_tree.Nnode < raxml_tree.Nnode

    meta = DataFrame({"taxa": sample_tree.tip_label, "grp": ["G1" if i < 6 else "G2" for i in range(13)]})
    p_op = ggtree(groupClade(groupOTU(sample_tree, ["A", "B", "C"]), mrca_de)) % meta
    p_op = p_op + geom_tiplab(aes(label="padded"), data=td_mutate(padded=lambda r: label_pad(r["label"], pad=".")))
    assert is_ggplot(p_op % rtree(13))

    # 6. ggdensitree, tanglegram / ggdoubletree, and multiplot / patchwork
    trees = read_tree(system_file("extdata", "ggdensitree_example.tree", package="ggtree"))
    outputs_to_render["11_ggdensitree"] = (
        ggdensitree(trees, alpha=0.3, colour="steelblue") + geom_tiplab(size=2.5) + hexpand(0.3) + ggtitle("Densitree"),
        (9.0, 6.5),
    )

    t1 = read_tree(text="(((a:1,b:1):1,c:2):1,(d:1.5,e:1.5):1.5);")
    t2 = read_tree(text="((e:1.5,d:1.5):1.5,(c:2,(b:1,a:1):1):1);")
    assoc = DataFrame({"left": ["a", "b", "c", "d", "e"], "right": ["a", "b", "c", "d", "e"]})
    tg = as_tanglegram(t1, t2, assoc=assoc)
    tg_df = fortify(tg, optimize=True, optimize_side="both")
    assert tg_df.attrs["tangle_optimize"]["after"] == 0
    outputs_to_render["12_ggdoubletree"] = (
        ggdoubletree(
            ggtree(t1) + geom_tiplab(),
            ggtree(t2) + geom_tiplab(),
            assoc=assoc,
            optimize=True,
            optimize_side="both",
            link_color="steelblue",
        )
        + ggtitle("Optimized Tanglegram (0 crossings)"),
        (9.5, 6.0),
    )

    p_multi = multiplot(
        ggtree(tr10) + geom_tiplab(size=2.5) + theme_tree() + ggtitle("Rectangular"),
        ggtree(tr10, layout="circular") + geom_tiplab2(size=2.5) + ggtitle("Circular"),
        ncol=2,
    )
    p_patch = (
        (ggtree(tr10) + geom_tiplab(size=2.5) + ggtitle("Rectangular"))
        | (ggtree(tr10, layout="fan", open_angle=120) + ggtitle("Fan (open=120)"))
    ) / (ggtree(tr10, layout="tidy") + geom_tiplab(size=2.5) + ggtitle("Tidy"))
    outputs_to_render["13_multiplot"] = (p_multi, (10.0, 5.5))
    outputs_to_render["14_patchwork"] = (p_patch, (10.0, 8.0))
    assert is_ggplot(plot_list(ggtree(tr10), ggtree(tr10, layout="circular")))

    # 7. Persist every composite figure to examples/output/ as SVG and PNG and validate headers
    out_dir = _REPO_ROOT / "examples" / "output"
    saved = _verify_and_save_outputs(outputs_to_render, out_dir)

    print(f"Saved and verified {len(saved)} figures (SVG + PNG) in {out_dir}:")
    for name, svg_p, png_p in saved:
        print(f"  {name:<30s} SVG={svg_p.stat().st_size:>7d}B  PNG={png_p.stat().st_size:>7d}B")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
