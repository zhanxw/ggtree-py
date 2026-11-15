"""
End-to-end showcase demonstrating True Vector PDF export and Interactive HTML /
Jupyter figures in ggtree-py (zero external runtime dependencies beyond numpy).

Run from the repository root with:
    python3 -s examples/interactive_pdf_example.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from ggtree import (
    DataFrame,
    aes,
    geom_tiplab,
    geom_tippoint,
    ggsave,
    ggtitle,
    ggtree,
    gheatmap,
    girafe,
    groupClade,
    hexpand,
    new_scale_fill,
    read_tree,
    scale_color_brewer,
    scale_fill_brewer,
    scale_fill_viridis_c,
    system_file,
    theme_tree2,
)
from ggtree.export import svg_to_pdf_pure_python


def main() -> int:
    tr = read_tree(system_file("extdata", "sample.nwk", package="treeio"))
    tips = list(tr.tip_label)  # 13 tips: A..M

    tr_clade = groupClade(tr, [17, 21])
    meta = DataFrame(
        {
            "label": tips,
            "habitat": [
                "Marine", "Marine", "Freshwater", "Terrestrial",
                "Terrestrial", "Freshwater", "Marine", "Freshwater",
                "Terrestrial", "Marine", "Terrestrial", "Freshwater", "Marine",
            ],
            "score": [12.4, 28.1, 9.5, 34.0, 18.2, 22.6, 41.3, 15.8, 26.0, 11.1, 31.5, 22.9, 16.4],
        }
    )

    p_tree = (
        ggtree(tr_clade, aes(color="group"))
        % meta
        + geom_tippoint(aes(fill="habitat", size="score"), shape=21, color="#1F2937")
        + geom_tiplab(size=3.0, offset=0.55)
        + scale_color_brewer(palette="Set1")
        + scale_fill_brewer(palette="Dark2")
        + hexpand(0.18)
        + theme_tree2()
        + ggtitle("Interactive Phylogeny & Trait Heatmap")
    )

    hm_mat = DataFrame(
        {
            "Temp": [round(10.0 + i * 1.4, 2) for i in range(len(tips))],
            "Salinity": [round(35.0 - (i % 5) * 6.2, 2) for i in range(len(tips))],
            "Depth": [round(5.0 + (i % 4) * 18.5, 2) for i in range(len(tips))],
        },
        rownames=tips,
    )
    p = (
        gheatmap(
            p_tree + new_scale_fill(),
            hm_mat,
            offset=6.5,
            width=0.35,
            colnames=True,
            colnames_angle=25,
            colnames_offset_y=-0.35,
            font_size=2.8,
        )
        + scale_fill_viridis_c(option="viridis")
        + hexpand(0.12)
    )

    out_dir = _REPO_ROOT / "examples" / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "16_interactive_and_pdf.pdf"
    html_path = out_dir / "16_interactive_and_pdf.html"

    # 1. True Vector PDF export (uses cairosvg -> rsvg-convert -> pure-Python vector PDF compiler)
    ggsave(pdf_path, plot=p, width=10.0, height=6.5)

    # Also exercise the pure-Python vector PDF fallback directly in memory
    pure_pdf = svg_to_pdf_pure_python(p.to_svg(width=1000, height=650), page_width_pt=720.0, page_height_pt=468.0)

    # 2. Interactive HTML / Jupyter widget export
    iv = girafe(p, width_svg=10.0, height_svg=6.5, dpi=100, title="ggtree-py Interactive Phylogeny")
    iv.save(html_path, width=10.0, height=6.5, dpi=100)

    pdf_bytes = pdf_path.read_bytes()
    html_txt = html_path.read_text(encoding="utf-8")
    iframe_repr = iv._repr_html_()

    assert pdf_bytes.startswith(b"%PDF-"), "Expected valid PDF header"
    assert pure_pdf.startswith(b"%PDF-1.4"), "Expected pure-Python %PDF-1.4 header"
    assert "ggtree-node-overlay" in html_txt
    assert "data-tooltip=" in html_txt
    assert iframe_repr.startswith('<iframe srcdoc="')

    print(
        f"Saved {pdf_path} ({len(pdf_bytes)} bytes), {html_path} ({len(html_txt)} bytes); "
        f"pure-Python PDF = {len(pure_pdf)} bytes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
