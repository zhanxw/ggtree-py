"""
End-to-end showcase demonstrating pandas integration with ggtree-py.

Run from the repository root with:
    python3 -s examples/pandas_example.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np

try:
    import pandas as pd
except ImportError as exc:
    raise SystemExit(
        "This example requires pandas (`pip install pandas` or `pip install 'ggtree-py[pandas]'`)."
    ) from exc

from ggtree import (
    aes,
    facet_labeller,
    facet_plot,
    facet_widths,
    fortify,
    geom_col,
    geom_tiplab,
    geom_tippoint,
    ggsave,
    ggtitle,
    ggtree,
    gheatmap,
    groupOTU,
    hexpand,
    inset,
    is_ggplot,
    nodepie,
    read_beast,
    read_tree,
    scale_color_brewer,
    scale_fill_brewer,
    scale_fill_viridis_c,
    scale_size,
    system_file,
    theme_tree2,
    xlim_expand,
    xlim_tree,
)


def main() -> int:
    nwk_file = system_file("extdata", "sample.nwk", package="treeio")
    sample_tree = read_tree(nwk_file)
    tips = list(sample_tree.tip_label)  # 13 tips: A..M

    # ------------------------------------------------------------------
    # Panel A: p % pdf + pd.Categorical + nullable Float64 / pd.NA +
    #          lambda d: d.to_pandas().query(...)
    # ------------------------------------------------------------------
    habitats = [
        "Marine",
        "Marine",
        "Freshwater",
        "Terrestrial",
        "Terrestrial",
        "Freshwater",
        "Marine",
        "Freshwater",
        "Terrestrial",
        "Marine",
        "Terrestrial",
        "Freshwater",
        "Marine",
    ]
    longevities = [12.0, 28.0, 9.0, 34.0, pd.NA, 19.0, 41.0, 15.0, 26.0, 11.0, 31.0, 22.0, 16.0]
    meta_pdf = pd.DataFrame(
        {
            "habitat": pd.Categorical(
                habitats,
                categories=["Marine", "Freshwater", "Terrestrial"],
                ordered=True,
            ),
            "longevity": pd.array(longevities, dtype="Float64"),
        },
        index=tips,
    )

    p1 = (
        ggtree(sample_tree)
        % meta_pdf
        + geom_tippoint(
            data=lambda d: d.to_pandas().query("isTip and longevity >= 25"),
            color="#F39C12",
            size=8.2,
            alpha=0.40,
        )
        + geom_tippoint(aes(color="habitat", size="longevity"), alpha=0.92)
        + geom_tiplab(aes(color="habitat"), offset=0.85, size=2.8)
        + scale_color_brewer(palette="Set1")
        + scale_size(range=(2.2, 5.5))
        + hexpand(0.32)
        + theme_tree2()
        + ggtitle("A. p % pd.DataFrame (index + pd.NA + query)")
    )

    # ------------------------------------------------------------------
    # Panel B: pd.DataFrame / pd.read_csv directly into gheatmap(...)
    # ------------------------------------------------------------------
    genotype_file = system_file("extdata/Genotype.txt", package="ggtree")
    genotype_pdf = pd.read_csv(genotype_file, sep="\t", index_col=0)
    beast_tree = read_beast(system_file("extdata/MCC_FluA_H3.tree", package="ggtree"))
    assert is_ggplot(gheatmap(ggtree(beast_tree, mrsd="2013-01-01"), genotype_pdf, offset=5, width=0.5))

    rng = np.random.default_rng(7)
    expr_pdf = pd.DataFrame(
        {
            "GeneA": np.round(rng.uniform(0.5, 4.5, size=len(tips)), 2),
            "GeneB": np.round(rng.uniform(0.2, 4.8, size=len(tips)), 2),
            "GeneC": np.round(rng.uniform(0.8, 5.0, size=len(tips)), 2),
            "GeneD": np.round(rng.uniform(0.1, 4.2, size=len(tips)), 2),
        },
        index=tips,
    )
    p2_base = ggtree(sample_tree) + geom_tiplab(size=2.6, offset=0.25)
    p2 = (
        gheatmap(
            p2_base,
            expr_pdf,
            offset=3.2,
            width=0.45,
            font_size=2.6,
            colnames_angle=-30,
            hjust=0.0,
        )
        + scale_fill_viridis_c(option="D")
        + ggtitle("B. gheatmap(p, expr_pdf)")
    )

    # ------------------------------------------------------------------
    # Panel C: facet_plot + nodepie directly from pandas.DataFrames
    # ------------------------------------------------------------------
    trait_pdf = pd.DataFrame(
        {
            "label": tips,
            "score": [3.2, 6.8, 2.5, 8.1, 4.4, 5.0, 9.2, 3.8, 7.0, 2.9, 7.6, 5.4, 4.1],
            "habitat": pd.Categorical(
                habitats,
                categories=["Marine", "Freshwater", "Terrestrial"],
            ),
        }
    )
    internal_nodes = list(range(len(tips) + 1, len(tips) + sample_tree.Nnode + 1))[:5]
    pie_pdf = pd.DataFrame(
        {
            "Marine": [0.7, 0.2, 0.5, 0.1, 0.6],
            "Freshwater": [0.2, 0.6, 0.3, 0.2, 0.2],
            "Terrestrial": [0.1, 0.2, 0.2, 0.7, 0.2],
        },
        index=internal_nodes,
    )
    pies = nodepie(
        pie_pdf,
        cols=["Marine", "Freshwater", "Terrestrial"],
        color={"Marine": "#E41A1C", "Freshwater": "#377EB8", "Terrestrial": "#4DAF4A"},
        alpha=0.85,
    )
    p3 = inset(ggtree(sample_tree) + geom_tiplab(size=2.6, offset=0.4), pies, width=0.065, height=0.065)
    p3 = facet_plot(
        p3,
        panel="TraitScore",
        data=trait_pdf,
        geom=geom_col,
        mapping=aes(x="score", fill="habitat"),
    )
    p3 = (
        facet_widths(
            facet_labeller(p3, {"Tree": "Tree + Pies", "TraitScore": "Trait Score"}),
            {"Tree": 0.62, "TraitScore": 0.38},
        )
        + xlim_tree(66)
        + xlim_expand((0, 10), "TraitScore")
        + scale_fill_brewer(palette="Set1")
        + ggtitle("C. facet_plot + nodepie(pd.DataFrame)")
    )

    # ------------------------------------------------------------------
    # Panel D: Round-trip fortify(tr).to_pandas() -> pandas ops -> ggtree(pdf)
    # ------------------------------------------------------------------
    tbl_pdf = fortify(sample_tree).to_pandas()
    focus_tips = pd.Series(["I", "J", "K", "L", "M"])
    grouped_tbl = groupOTU(tbl_pdf, focus_tips, group_name="clade_focus").to_pandas(categorical=True)
    grouped_tbl["depth_band"] = np.where(
        grouped_tbl["x"] >= grouped_tbl["x"].median(),
        "Distal",
        "Basal",
    )
    p4 = (
        ggtree(grouped_tbl, mapping=aes(color="clade_focus"))
        + geom_tippoint(
            grouped_tbl[grouped_tbl["isTip"] & (grouped_tbl["depth_band"] == "Distal")],
            aes(x="x", y="y"),
            color="#8E44AD",
            size=3.6,
            inherit_aes=False,
        )
        + geom_tiplab(size=2.8, offset=0.55)
        + scale_color_brewer(palette="Dark2")
        + hexpand(0.28)
        + theme_tree2()
        + ggtitle("D. fortify().to_pandas() -> ggtree(pdf)")
    )

    composite = (p1 | p2) / (p3 | p4)
    out_dir = _REPO_ROOT / "examples" / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    svg_path = out_dir / "15_pandas_integration.svg"
    png_path = out_dir / "15_pandas_integration.png"
    ggsave(svg_path, plot=composite, width=14.0, height=10.5, dpi=120)
    ggsave(png_path, plot=composite, width=14.0, height=10.5, dpi=120)

    svg_txt = svg_path.read_text(encoding="utf-8")
    png_bytes = png_path.read_bytes()
    assert svg_txt.lstrip().startswith("<svg")
    assert png_bytes.startswith(b"\x89PNG\r\n\x1a\n")
    print(f"Saved {svg_path} ({len(svg_txt)} bytes) and {png_path} ({len(png_bytes)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
