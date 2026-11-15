"""
Tests for True PDF vector export (including pure-Python fallback) and
interactive HTML / Jupyter figure rendering (to_html, to_interactive, girafe).
"""

from __future__ import annotations

import os
import re
import tempfile
import unittest
import zlib
from pathlib import Path

from ggtree import (
    DataFrame,
    InteractivePlot,
    aes,
    facet_plot,
    geom_col,
    geom_tiplab,
    geom_tippoint,
    ggsave,
    ggplot2,
    ggtree,
    gheatmap,
    girafe,
    groupClade,
    hexpand,
    read_tree,
    rtree,
    scale_fill_viridis_c,
    svg_to_pdf_bytes,
    system_file,
    theme_tree2,
)
from ggtree.export import svg_to_pdf_pure_python


class TestExportInteractiveAndPDF(unittest.TestCase):
    def setUp(self):
        self.tr = read_tree(text="((t1:1.0,t2:1.5)n1:0.8,(t3:1.2,(t4:0.5,t5:0.7)n3:0.9)n2:1.1)root;")
        self.meta = DataFrame(
            {
                "node": [1, 2, 3, 4, 5],
                "habitat": ["Marine", "Freshwater", "Marine", "Terrestrial", "Freshwater"],
                "score": [10.5, 22.0, 15.25, 31.0, 8.5],
            }
        )

    @staticmethod
    def _extract_pure_pdf_stream(pdf_bytes: bytes) -> str:
        m = re.search(rb"stream\n(.*?)\nendstream", pdf_bytes, re.S)
        if not m:
            return ""
        return zlib.decompress(m.group(1)).decode("latin1")

    def test_to_pdf_rectangular_and_polar(self):
        p_rect = ggtree(self.tr) % self.meta + geom_tippoint(aes(color="habitat")) + geom_tiplab(offset=0.1) + hexpand(0.2)
        pdf_rect = p_rect.to_pdf(width=800, height=600)
        self.assertIsInstance(pdf_rect, bytes)
        self.assertTrue(pdf_rect.startswith(b"%PDF-"))
        self.assertIn(b"%%EOF", pdf_rect[-32:])
        self.assertGreater(len(pdf_rect), 1000)

        p_circ = ggtree(self.tr, layout="circular") + geom_tippoint(color="#2563EB") + geom_tiplab()
        pdf_circ = p_circ.to_pdf(width=700, height=700)
        self.assertTrue(pdf_circ.startswith(b"%PDF-"))
        self.assertGreater(len(pdf_circ), 1000)

        pdf_pt = p_rect.to_pdf(width=800, height=600, page_width_pt=648.0, page_height_pt=432.0)
        self.assertTrue(pdf_pt.startswith(b"%PDF-"))

    def test_pure_python_svg_to_pdf_fallback(self):
        p_rect = ggtree(self.tr) + geom_tippoint() + geom_tiplab()
        p_circ = ggtree(self.tr, layout="circular") + geom_tiplab()
        mat = DataFrame(
            {"TraitA": [1.0, 2.5, 0.5, 3.0, 1.8], "TraitB": [2.0, 1.0, 3.2, 0.4, 2.7]},
            rownames=list(self.tr.tip_label),
        )
        p_hm = gheatmap(ggtree(self.tr) + geom_tiplab(offset=0.05), mat, offset=0.8, width=0.5)
        pl = p_rect | p_circ

        for plot_obj, desc in (
            (p_rect, "rectangular"),
            (p_circ, "circular"),
            (p_hm, "gheatmap"),
            (pl, "plotlist"),
        ):
            svg = plot_obj.to_svg(width=760, height=540)
            pdf_b = svg_to_pdf_pure_python(svg, page_width_pt=540.0, page_height_pt=380.0)
            self.assertTrue(pdf_b.startswith(b"%PDF-1.4"), msg=f"{desc} header mismatch")
            self.assertIn(b"/MediaBox [0 0 540.00 380.00]", pdf_b, msg=f"{desc} MediaBox mismatch")
            self.assertIn(b"/Filter /FlateDecode", pdf_b, msg=f"{desc} FlateDecode missing")
            self.assertTrue(pdf_b.endswith(b"%%EOF\n"), msg=f"{desc} EOF missing")
            stream_txt = self._extract_pure_pdf_stream(pdf_b)
            self.assertIn("BT", stream_txt, msg=f"{desc} text operator missing")
            self.assertIn("Tj", stream_txt, msg=f"{desc} Tj operator missing")
            self.assertIn("t1", stream_txt, msg=f"{desc} tip label missing from PDF stream")
            if desc == "circular":
                self.assertIn(" l\n", stream_txt, "circular arc polyline segments missing")
                self.assertIn(" cm\n", stream_txt, "rotated circular tip label matrix missing")

        # Directly verify SVG arc ('A'), cubic/quadratic curves ('C','Q'), and opacity in pure-Python PDF
        bezier_svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="150" viewBox="0 0 200 150">'
            '<path d="M 20 80 A 40 40 0 0 1 100 80 Q 140 20 180 80" fill="none" stroke="#2563EB" stroke-width="2" stroke-opacity="0.7"/>'
            '</svg>'
        )
        bezier_pdf = svg_to_pdf_pure_python(bezier_svg, page_width_pt=200.0, page_height_pt=150.0)
        bezier_stream = self._extract_pure_pdf_stream(bezier_pdf)
        self.assertIn(" c\n", bezier_stream, "SVG arc/curve cubic Bezier missing")
        self.assertIn("/GS1 gs", bezier_pdf.decode("latin1") + bezier_stream)

        # Also test svg_to_pdf_bytes wrapper
        wrapped = svg_to_pdf_bytes(p_rect.to_svg(width=640, height=480), width=640, height=480)
        self.assertTrue(wrapped.startswith(b"%PDF-"))

    def test_plotlist_to_pdf_and_to_interactive(self):
        p1 = ggtree(self.tr) + geom_tiplab()
        p2 = ggtree(self.tr, layout="slanted") + geom_tippoint(color="#16A34A") + geom_tiplab()
        pl = p1 | p2

        pdf_b = pl.to_pdf(width=900, height=450)
        self.assertTrue(pdf_b.startswith(b"%PDF-"))

        iv = pl.to_interactive(width=900, height=450, title="Side-by-Side Trees")
        self.assertIsInstance(iv, InteractivePlot)
        html_str = pl.to_html(width=900, height=450, interactive=True, title="Side-by-Side Trees")
        self.assertIn("Side-by-Side Trees", html_str)
        self.assertGreaterEqual(html_str.count('class="ggtree-node-overlay"'), 2)

    def test_ggsave_formats_device_and_validation(self):
        p = ggtree(self.tr) + geom_tiplab() + geom_tippoint()
        with tempfile.TemporaryDirectory() as tmpdir:
            td = Path(tmpdir)

            # Explicit extensions
            svg_p = td / "tree.svg"
            png_p = td / "tree.png"
            pdf_p = td / "tree.pdf"
            html_iv_p = td / "tree_iv.html"
            html_st_p = td / "tree_static.html"

            ggsave(svg_p, plot=p, width=7.5, height=5.0, dpi=100)
            ggsave(png_p, plot=p, width=7.5, height=5.0, dpi=100)
            ggsave(pdf_p, plot=p, width=7.5, height=5.0)
            ggsave(html_iv_p, plot=p, width=8.0, height=5.5, interactive=True)
            ggsave(html_st_p, plot=p, width=8.0, height=5.5, interactive=False)

            self.assertTrue(svg_p.read_text(encoding="utf-8").lstrip().startswith("<svg"))
            self.assertNotIn("ggtree-node-overlay", svg_p.read_text(encoding="utf-8"))
            self.assertTrue(png_p.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertTrue(pdf_p.read_bytes().startswith(b"%PDF-"))

            iv_txt = html_iv_p.read_text(encoding="utf-8")
            self.assertIn("ggtree-interactive-widget", iv_txt)
            self.assertIn("ggtree-node-overlay", iv_txt)

            st_txt = html_st_p.read_text(encoding="utf-8")
            self.assertTrue(st_txt.startswith("<!DOCTYPE html>"))
            self.assertNotIn("ggtree-node-overlay", st_txt)
            self.assertNotIn("ggtree-interactive-widget", st_txt)

            # device="pdf" overrides extensionless or custom suffix
            dev_p = td / "custom_output.dat"
            ggsave(dev_p, plot=p, width=6.0, height=4.0, device="pdf")
            self.assertTrue(dev_p.read_bytes().startswith(b"%PDF-"))

            # Unsupported extension raises ValueError
            with self.assertRaises(ValueError):
                ggsave(td / "bad.xyz", plot=p)

            # Uses last_plot() when plot=None
            last_p = td / "last_plot.pdf"
            ggsave(last_p, width=6.0, height=4.0)
            self.assertTrue(last_p.read_bytes().startswith(b"%PDF-"))

    def test_to_html_interactive_overlay_vs_clean_to_svg(self):
        p = (
            ggtree(self.tr) % self.meta
            + geom_tippoint(aes(color="habitat", size="score"))
            + geom_tiplab(offset=0.1)
            + theme_tree2()
        )
        clean_svg = p.to_svg(width=800, height=560)
        self.assertNotIn("ggtree-node-overlay", clean_svg)
        self.assertNotIn("ggtree-hit-node", clean_svg)
        self.assertNotIn("data-tooltip=", clean_svg)

        html_iv = p.to_html(width=820, height=580, interactive=True, title="Annotated Tree")
        for token in (
            "Annotated Tree",
            "ggtree-interactive-widget",
            "ggtree-node-overlay",
            "ggtree-hit-node",
            "ggtree-clade-edge",
            "data-meta=",
            "Terrestrial",
            "Marine",
            "toggleCollapse",
            "selectClade",
            "applySearch",
            "dl-svg",
            "dl-png",
        ):
            self.assertIn(token, html_iv)

        # Static HTML & fragment HTML
        html_static = p.to_html(width=800, height=560, interactive=False)
        self.assertTrue(html_static.startswith("<!DOCTYPE html>"))
        self.assertIn("<svg", html_static)
        self.assertNotIn("ggtree-node-overlay", html_static)

        frag = p.to_html(width=800, height=560, interactive=True, full_html=False)
        self.assertFalse(frag.lstrip().startswith("<!DOCTYPE html>"))
        self.assertIn("<style>", frag)
        self.assertIn("ggtree-interactive-widget", frag)
        self.assertIn("<script>", frag)

    def test_gheatmap_and_facet_plot_tooltips(self):
        mat = DataFrame(
            {
                "Gene1": [0.1, 2.4, 1.1, 3.8, 0.6],
                "Gene2": [3.2, 0.9, 2.7, 1.5, 4.1],
            },
            rownames=list(self.tr.tip_label),
        )
        p_hm = gheatmap(ggtree(self.tr) + geom_tiplab(offset=0.05), mat, offset=0.8, width=0.5) + scale_fill_viridis_c()
        clean_hm_svg = p_hm.to_svg(width=840, height=560)
        self.assertNotIn("data-tooltip=", clean_hm_svg)
        self.assertNotIn("ggtree-node-overlay", clean_hm_svg)

        hm_html = p_hm.to_html(width=840, height=560, interactive=True)
        self.assertIn("data-tooltip=", hm_html)
        self.assertIn("col=Gene1", hm_html)
        self.assertIn("value=", hm_html)
        self.assertIn("ggtree-node-overlay", hm_html)

        dd = DataFrame(
            {
                "id": list(self.tr.tip_label),
                "abundance": [12.5, 28.0, 17.5, 34.0, 9.0],
            }
        )
        p_fp = facet_plot(ggtree(self.tr) + geom_tiplab(), panel="Abundance", data=dd, geom=geom_col, mapping=aes(x="abundance"))
        clean_fp_svg = p_fp.to_svg(width=840, height=560)
        self.assertNotIn("data-tooltip=", clean_fp_svg)

        fp_html = p_fp.to_html(width=840, height=560, interactive=True)
        self.assertIn("data-tooltip=", fp_html)
        self.assertIn("panel=Abundance", fp_html)
        self.assertIn("x=12.5", fp_html)

    def test_interactive_plot_girafe_repr_and_save(self):
        p = ggtree(self.tr) + geom_tiplab() + geom_tippoint(color="#D97706")
        iv = p.to_interactive(width=820, height=580, title="Jupyter Widget")
        self.assertIsInstance(iv, InteractivePlot)

        iv_g = girafe(p, width_svg=8.2, height_svg=5.8, dpi=100, title="Girafe Widget")
        self.assertIsInstance(iv_g, InteractivePlot)
        self.assertEqual(iv_g.width, 820)
        self.assertEqual(iv_g.height, 580)
        self.assertIs(ggplot2.girafe(p).plot, p)

        iframe_html = iv._repr_html_()
        self.assertTrue(iframe_html.startswith('<iframe srcdoc="'))
        self.assertIn('sandbox="allow-scripts allow-downloads"', iframe_html)
        self.assertIn("Jupyter Widget", iframe_html)

        self.assertNotIn("ggtree-node-overlay", iv.to_svg())
        self.assertTrue(iv.to_png().startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertTrue(iv.to_pdf().startswith(b"%PDF-"))

        with tempfile.TemporaryDirectory() as tmpdir:
            td = Path(tmpdir)
            out_html = td / "widget.html"
            out_pdf = td / "widget.pdf"
            iv.save(out_html, width=8.2, height=5.8, dpi=100)
            iv.save(out_pdf, width=8.2, height=5.8)
            self.assertIn("ggtree-interactive-widget", out_html.read_text(encoding="utf-8"))
            self.assertTrue(out_pdf.read_bytes().startswith(b"%PDF-"))

    def test_sample_nwk_interactive_and_pdf(self):
        tr = read_tree(system_file("extdata", "sample.nwk", package="treeio"))
        tr_grp = groupClade(tr, [15, 19])
        p = ggtree(tr_grp, aes(color="group")) + geom_tippoint(size=3.0) + geom_tiplab(offset=0.2) + hexpand(0.2)
        self.assertTrue(p.to_pdf(width=760, height=520).startswith(b"%PDF-"))
        html_doc = p.to_html(width=760, height=520, interactive=True)
        self.assertIn("group", html_doc)
        self.assertIn("ggtree-hit-node", html_doc)


if __name__ == "__main__":
    unittest.main()
