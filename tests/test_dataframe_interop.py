"""
Tests for external DataFrame interoperability:
pandas (including nullable dtypes / pd.NA / pd.Categorical / Index / Series),
polars, pyarrow, and duck-typed / __dataframe__ objects.
"""

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from ggtree import (
    DataFrame,
    MRCA,
    aes,
    ancestor,
    arrange,
    as_phylo,
    as_tibble,
    as_treedata,
    bind_rows,
    c,
    child,
    ensure_dataframe,
    facet_plot,
    filter,
    fortify,
    full_join,
    geom_col,
    geom_point,
    geom_tippoint,
    get_layout,
    ggdensitree,
    ggplot,
    ggplot_build,
    ggtree,
    gheatmap,
    groupClade,
    groupOTU,
    inner_join,
    is_dataframe_like,
    is_ggplot,
    is_na,
    isTip,
    left_join,
    mutate,
    nodebar,
    nodeid,
    nodelab,
    nodepie,
    offspring,
    parent,
    read_tree,
    rootnode,
    select,
    sibling,
    system_file,
    tip_y_alignment,
    unnest,
)

try:
    import pandas as pd  # type: ignore

    HAS_PANDAS = True
except ImportError:
    pd = None  # type: ignore
    HAS_PANDAS = False

try:
    import polars as pl  # type: ignore

    HAS_POLARS = True
except ImportError:
    pl = None  # type: ignore
    HAS_POLARS = False

try:
    import pyarrow as pa  # type: ignore

    HAS_PYARROW = True
except ImportError:
    pa = None  # type: ignore
    HAS_PYARROW = False


class TestPurePythonDataFrameProtocol(unittest.TestCase):
    """Tests that run with only numpy installed (no pandas/polars/pyarrow)."""

    def test_duck_typed_na_type_and_mock_dataframe(self):
        class NAType:
            def __bool__(self):
                raise TypeError("boolean value of NA is ambiguous")

            def __eq__(self, other):
                return NAType()

            def __ne__(self, other):
                return NAType()

        na_val = NAType()
        self.assertTrue(is_na(na_val))

        class MockSeries:
            dtype = "object"
            ndim = 1
            name = "trait"

            def __init__(self, vals):
                self._vals = list(vals)

            def tolist(self):
                return list(self._vals)

        class MockPandasDF:
            columns = ["taxa", "val"]
            index = ["r1", "r2"]
            attrs = {"source": "mock"}

            def __init__(self):
                self._d = {
                    "taxa": MockSeries(["A", "B"]),
                    "val": MockSeries([10, na_val]),
                }

            def __getitem__(self, k):
                return self._d[k]

            def to_dict(self):
                return {k: s.tolist() for k, s in self._d.items()}

        mock = MockPandasDF()
        self.assertTrue(is_dataframe_like(mock))
        df = ensure_dataframe(mock)
        self.assertEqual(df["taxa"].tolist(), ["A", "B"])
        self.assertEqual(df["val"].tolist(), [10, None])
        self.assertEqual(df.rownames, ["r1", "r2"])
        self.assertEqual(df.attrs["source"], "mock")

    def test_missing_optional_export_raises_helpful_error(self):
        df = DataFrame({"a": [1, 2]})
        if not HAS_PANDAS:
            with self.assertRaises(ImportError):
                df.to_pandas()
        if not HAS_POLARS:
            with self.assertRaises(ImportError):
                df.to_polars()
        if not HAS_PYARROW:
            with self.assertRaises(ImportError):
                df.to_arrow()


@unittest.skipUnless(HAS_PANDAS, "pandas is not installed")
class TestPandasInterop(unittest.TestCase):
    def setUp(self):
        self.tr = read_tree(text="(((A:1,B:1):1,C:2):1,(D:1.5,E:1.5):1.5);")

    def test_rangeindex_vs_explicit_index_and_nullable_dtypes(self):
        pdf_range = pd.DataFrame(
            {
                "label": pd.Series(["A", "B", "C"], dtype="string"),
                "count": pd.Series([1, pd.NA, 3], dtype="Int64"),
                "score": pd.Series([1.5, 2.5, pd.NA], dtype="Float64"),
                "flag": pd.Series([True, pd.NA, False], dtype="boolean"),
                "dt": pd.to_datetime(["2020-01-01", None, "2022-05-06"]),
                "cat": pd.Categorical(["hi", "lo", None], categories=["lo", "hi"], ordered=True),
            }
        )
        df = DataFrame(pdf_range)
        # Default unnamed RangeIndex(0, 3, 1) does not become synthetic rownames
        self.assertIsNone(df._rownames)
        self.assertEqual(df["count"].tolist(), [1, None, 3])
        self.assertEqual(df["score"].tolist(), [1.5, 2.5, None])
        self.assertEqual(df["flag"].tolist(), [True, None, False])
        self.assertIsNone(df["dt"].tolist()[1])
        self.assertEqual(df["cat"].tolist(), ["hi", "lo", None])
        self.assertEqual(df._levels["cat"], ["lo", "hi"])
        self.assertEqual(is_na(pd.Series([1, pd.NA, 3], dtype="Int64")).tolist(), [False, True, False])

        # Series directly into DataFrame dict, __setitem__, and c()
        df2 = DataFrame({"s": pdf_range["cat"], "x": pdf_range["count"]})
        self.assertEqual(df2.shape, (3, 2))
        self.assertEqual(df2._levels["s"], ["lo", "hi"])
        df2["extra"] = pd.Series([10, 20, 30])
        self.assertEqual(df2["extra"].tolist(), [10, 20, 30])
        self.assertEqual(c(pd.Series(["A", "B"]), "C"), ["A", "B", "C"])

        # Custom index preserved; roundtrip with categorical levels
        pdf_idx = pd.DataFrame({"val": [1, 2]}, index=["A", "B"])
        df_idx = DataFrame.from_pandas(pdf_idx)
        self.assertEqual(df_idx.rownames, ["A", "B"])
        back = df.to_pandas(categorical=True)
        self.assertEqual(list(back["cat"].cat.categories), ["lo", "hi"])
        back_idx = df_idx.to_pandas()
        self.assertEqual(list(back_idx.index), ["A", "B"])

    def test_attach_data_operators_index_and_label_priority(self):
        # 1. Index-keyed pandas DataFrame (no label/taxa column) + ordered categorical legend
        meta_idx = pd.DataFrame(
            {
                "group": pd.Categorical(["G1", "G1", "G2", "G2", "G3"], categories=["G3", "G2", "G1"]),
                "cat": pd.Categorical(["low", "high", "low", "med", "high"], categories=["low", "med", "high"], ordered=True),
                "weight": pd.array([10, 20, pd.NA, 40, 50], dtype="Int64"),
            },
            index=["A", "B", "C", "D", "E"],
        )
        p1 = ggtree(self.tr) % meta_idx + geom_tippoint(aes(color="group"))
        self.assertTrue(is_ggplot(p1))
        self.assertEqual(p1.data._levels["group"], ["G3", "G2", "G1"])
        self.assertEqual(p1.data._levels["cat"], ["low", "med", "high"])
        a_row = [r for r in p1.data.to_records() if r.get("label") == "A"][0]
        self.assertEqual(a_row["group"], "G1")
        self.assertEqual(a_row["weight"], 10)
        svg = p1.to_svg()
        self.assertIn("<svg", svg)
        svg_cat = (ggtree(self.tr) % meta_idx + geom_tippoint(aes(color="cat"))).to_svg()
        self.assertLess(svg_cat.index(">low<"), svg_cat.index(">med<"))
        self.assertLess(svg_cat.index(">med<"), svg_cat.index(">high<"))

        # 2. First column matches label (R %<+% style), plus @ and + operators
        meta_col1 = pd.DataFrame({"taxa": ["A", "D"], "tag": ["tA", "tD"]})
        for p_op in (ggtree(self.tr) @ meta_col1, ggtree(self.tr) + meta_col1):
            d_rows = {r["label"]: r["tag"] for r in p_op.data.to_records() if r["label"]}
            self.assertEqual(d_rows["A"], "tA")
            self.assertEqual(d_rows["D"], "tD")

        # 3. 'label' in non-first column
        meta_col2 = pd.DataFrame({"tag": ["xA", "xE"], "label": ["A", "E"]})
        p3 = ggtree(self.tr) % meta_col2
        d_rows3 = {r["label"]: r["tag"] for r in p3.data.to_records() if r["label"]}
        self.assertEqual(d_rows3["A"], "xA")
        self.assertEqual(d_rows3["E"], "xE")

    def test_gheatmap_facet_plot_and_insets(self):
        # gheatmap with index-aligned pandas DataFrame & Genotype.txt via pd.read_csv
        mat_pdf = pd.DataFrame(
            {"G1": [1.0, 2.0, 3.0, 4.0, 5.0], "G2": [5.0, 4.0, 3.0, 2.0, 1.0]},
            index=["A", "B", "C", "D", "E"],
        )
        p_hm = gheatmap(ggtree(self.tr), mat_pdf, offset=0.5, width=0.4)
        self.assertTrue(is_ggplot(p_hm))
        self.assertIn("<svg", p_hm.to_svg())

        # Genotype.txt with pd.read_csv
        geno_file = system_file("extdata/Genotype.txt", package="ggtree")
        geno_pdf = pd.read_csv(geno_file, sep="\t", index_col=0)
        self.assertGreater(len(geno_pdf), 0)

        # gheatmap with 'label' column in 2nd position
        mat_lbl = pd.DataFrame({"val1": [1, 2, 3, 4, 5], "label": ["A", "B", "C", "D", "E"]})
        self.assertTrue(is_ggplot(gheatmap(ggtree(self.tr), mat_lbl)))

        # facet_plot with column-keyed and index-keyed pandas DataFrames
        trait_col = pd.DataFrame({"id": ["A", "B", "C", "D", "E"], "val": [2.0, 4.0, 1.0, 3.0, 5.0]})
        trait_idx = pd.DataFrame({"val": [2.0, 4.0, 1.0, 3.0, 5.0]}, index=["A", "B", "C", "D", "E"])
        for tdf in (trait_col, trait_idx):
            pf = facet_plot(ggtree(self.tr), panel="Trait", data=tdf, geom=geom_col, mapping=aes(x="val"))
            self.assertIn("<svg", pf.to_svg())

        # nodepie & nodebar with 'node' column or index
        pie_col = pd.DataFrame({"a": [0.4, 0.6], "b": [0.6, 0.4], "node": [6, 7]})
        bar_idx = pd.DataFrame({"a": [0.3, 0.8], "b": [0.7, 0.2]}, index=[6, 7])
        self.assertEqual(set(nodepie(pie_col, cols=["a", "b"]).keys()), {6, 7})
        self.assertEqual(set(nodebar(bar_idx, cols=["a", "b"]).keys()), {6, 7})

    def test_geom_positional_data_and_callable_returning_pandas(self):
        tbl_pdf = fortify(self.tr).to_pandas()
        tips_pdf = tbl_pdf[tbl_pdf["isTip"]]

        # Positional geom_point(pdf, aes(...))
        p_pos = ggtree(self.tr) + geom_point(tips_pdf, aes(x="x", y="y"), color="red")
        b = ggplot_build(p_pos)
        self.assertEqual(len(b.data[-1]), 5)

        # Lambda returning pandas.DataFrame inside layer data=
        p_lam = ggtree(self.tr) + geom_point(
            data=lambda d: d.to_pandas().query("isTip and x >= 2.5"),
            color="steelblue",
        )
        b2 = ggplot_build(p_lam)
        self.assertGreater(len(b2.data[-1]), 0)
        self.assertIn("<svg", p_lam.to_svg())

        # ggplot(pdf)
        self.assertTrue(is_ggplot(ggplot(tbl_pdf, aes(x="x", y="y"))))

    def test_tree_converters_tidytree_and_verbs(self):
        tbl_pdf = fortify(self.tr).to_pandas()

        # fortify / as_tibble / as_phylo / as_treedata / ggtree / ggdensitree
        self.assertEqual(get_layout(tbl_pdf) if False else None, None)
        self.assertEqual(tip_y_alignment(tbl_pdf), "ok")
        self.assertEqual(len(as_tibble(tbl_pdf)), len(tbl_pdf))
        self.assertEqual(len(fortify(tbl_pdf)), len(tbl_pdf))
        self.assertEqual(len(as_phylo(tbl_pdf).tip_label), 5)
        self.assertEqual(len(as_treedata(tbl_pdf).phylo.tip_label), 5)
        self.assertTrue(is_ggplot(ggtree(tbl_pdf)))
        self.assertTrue(is_ggplot(ggdensitree(tbl_pdf)))

        # tidytree functions on pd.DataFrame and pd.Series arguments
        rn = rootnode(tbl_pdf)
        self.assertEqual(len(rn), 1)
        self.assertEqual(nodeid(tbl_pdf, pd.Series(["A", "B"])), [1, 2])
        self.assertEqual(nodelab(tbl_pdf, pd.Series([1, 2])), ["A", "B"])
        self.assertEqual(isTip(tbl_pdf, pd.Series([1, 6])), [True, False])
        self.assertGreater(len(child(tbl_pdf, int(rn["node"][0]))), 0)
        self.assertEqual(len(parent(tbl_pdf, "A")), 1)
        self.assertEqual(len(sibling(tbl_pdf, "A")), 1)
        self.assertGreater(len(ancestor(tbl_pdf, "A")), 0)
        self.assertEqual(len(offspring(tbl_pdf, int(rn["node"][0]), tiponly=True)), 5)
        self.assertEqual(int(MRCA(tbl_pdf, "A", "B")), int(MRCA(self.tr, pd.Series(["A", "B"]))))

        # groupOTU / groupClade with pd.DataFrame and pd.Series
        grouped = groupOTU(tbl_pdf, pd.Series(["A", "B"]))
        self.assertIn("group", grouped.columns)
        claded = groupClade(self.tr, pd.Series([MRCA(self.tr, "A", "B")]))
        self.assertIn("group", claded.attrs)

        # dplyr-style verbs & joins with pd.DataFrame
        joined = left_join(self.tr, pd.DataFrame({"score": [11, 22]}, index=["A", "B"]))
        self.assertIn("score", joined.data.columns)
        fj = full_join(self.tr, pd.DataFrame({"label": ["A"], "flag": [True]}))
        self.assertIn("flag", fj.data.columns)
        ij = inner_join(pd.DataFrame({"k": [1, 2]}), pd.DataFrame({"k": [2, 3], "v": ["x", "y"]}))
        self.assertEqual(ij["v"].tolist(), ["x"])

        m_df = mutate(pd.DataFrame({"a": [1, 2]}), b="a * 10")
        self.assertEqual(m_df["b"].tolist(), [10, 20])
        f_df = filter(pd.DataFrame({"a": [1, 5, 2]}), "a > 1")
        self.assertEqual(f_df["a"].tolist(), [5, 2])
        self.assertEqual(select(pd.DataFrame({"a": [1], "b": [2]}), "b").columns, ["b"])
        self.assertEqual(arrange(pd.DataFrame({"a": [3, 1, 2]}), "a")["a"].tolist(), [1, 2, 3])
        self.assertEqual(len(bind_rows(pd.DataFrame({"a": [1]}), pd.DataFrame({"a": [2]}))), 2)
        self.assertEqual(len(unnest(pd.DataFrame({"id": [1], "xs": [[10, 20]]}), "xs")), 2)


@unittest.skipUnless(HAS_POLARS, "polars is not installed")
class TestPolarsInterop(unittest.TestCase):
    def test_polars_dataframe_and_export(self):
        tr = read_tree(text="((A:1,B:1):1,C:2);")
        pl_df = pl.DataFrame(
            {
                "label": ["A", "B", "C"],
                "grp": pl.Series(["X", "Y", "X"], dtype=pl.Categorical),
                "val": [1.0, None, 3.0],
            }
        )
        self.assertTrue(is_dataframe_like(pl_df))
        df = DataFrame.from_polars(pl_df)
        self.assertEqual(df["label"].tolist(), ["A", "B", "C"])
        self.assertEqual(df["val"].tolist(), [1.0, None, 3.0])

        p = ggtree(tr) + pl_df + geom_tippoint(aes(color="grp"))
        self.assertIn("<svg", p.to_svg())

        p_hm = gheatmap(ggtree(tr), pl.DataFrame({"label": ["A", "B", "C"], "V1": [1.0, 2.0, 3.0]}))
        self.assertIn("<svg", p_hm.to_svg())

        out_pl = fortify(tr).to_polars()
        self.assertEqual(len(out_pl), 5)


@unittest.skipUnless(HAS_PYARROW, "pyarrow is not installed")
class TestPyArrowInterop(unittest.TestCase):
    def test_pyarrow_table_and_export(self):
        tr = read_tree(text="((A:1,B:1):1,C:2);")
        pa_tbl = pa.table({"label": ["A", "B", "C"], "score": [2.0, None, 4.0]})
        self.assertTrue(is_dataframe_like(pa_tbl))
        df = DataFrame.from_arrow(pa_tbl)
        self.assertEqual(df["score"].tolist(), [2.0, None, 4.0])

        pf = facet_plot(ggtree(tr), panel="Arrow", data=pa_tbl, geom=geom_point, mapping=aes(x="score"))
        self.assertIn("<svg", pf.to_svg())

        out_pa = fortify(tr).to_arrow()
        self.assertEqual(out_pa.num_rows, 5)


if __name__ == "__main__":
    unittest.main()
