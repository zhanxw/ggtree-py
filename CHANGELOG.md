# Changelog

## 0.1.3 (2026-11-15)

- **Online Book website**: Linked the 13-chapter *treedata-book* Python edition ([`https://zhanxw.github.io/ggtree-py-book/`](https://zhanxw.github.io/ggtree-py-book/)) and 111-figure searchable gallery ([`https://zhanxw.github.io/ggtree-py-book/gallery.html`](https://zhanxw.github.io/ggtree-py-book/gallery.html)) in `pyproject.toml` (`Book` URL), `README.md`, and `docs/`.
- **DataFrame interoperability (`pandas`, `polars`, `pyarrow`, `__dataframe__`)**: Native duck-typed input support across `%` / `<<` / `@` / `+`, `gheatmap`, `facet_plot` / `geom_facet`, `nodepie` / `nodebar`, `geom_*` layers, `ggtree()`, `tidytree` helpers, and table verbs, plus `.to_pandas()`, `.to_polars()`, and `.to_arrow()` export methods (`examples/pandas_example.py`, `tests/test_dataframe_interop.py`).
- **Interactive HTML & Jupyter figures (`ggtree.export`)**: Added `to_html()`, `to_interactive()`, `girafe()`, and `InteractivePlot` (`_repr_html_()` `<iframe srcdoc=...>`) with pan/zoom, live search, node & heatmap/facet tile hover tooltips, click-to-highlight clades, double-click internal-node collapse toggle, and SVG/PNG download buttons.
- **True PDF vector export**: Added `to_pdf()`, `svg_to_pdf_bytes()`, and zero-dependency pure-Python `%PDF-1.4` vector compiler `svg_to_pdf_pure_python()`, and updated `ggsave()` to support `.pdf`, `.html`, `.svg`, `.png`, `device=`, and `interactive=` (`examples/interactive_pdf_example.py`, `tests/test_export_interactive_pdf.py`).

## 0.1.2 (2026-11-14)

- **Visual documentation gallery (28 figures)**: Added `docs/figures/*.png`, `docs/gallery.md`, and reproducible `docs/generate_doc_figures.py` (`120 DPI`), embedding rendered figures inline across `README.md`, `docs/index.md`, `docs/getting_started.md`, `docs/layouts_and_manipulation.md`, `docs/annotations_and_compositions.md`, and `docs/r_to_python_migration.md`.
- **SVG renderer & geom enhancements**: Added `geom_col`/`geom_bar` facet bar rendering, stacked and dodged `nodebar` insets, discrete/continuous mapped `aes(color=...)` and rescaled `aes(size=...)` (`scale_size`) support with deduplicated right-side facet legends, per-clade multi-fill support in `geom_balance`, `tiplab` deduplication in `ggdensitree`, and `node_label` preservation in `extract_clade()`.

## 0.1.1 (2026-11-14)

- **Hosted documentation & links**: Created public GitHub repository (`https://github.com/zhanxw/ggtree.py`) and GitHub Pages documentation site (`https://zhanxw.github.io/ggtree.py/`).
- **PyPI metadata & README links**: Replaced relative `docs/*.md`, `LICENSE`, and `CONDUCT.md` links with absolute GitHub Pages and GitHub URLs, and added `Homepage`, `Documentation`, `Repository`, `Issues`, `Changelog`, and `Upstream` URLs to `pyproject.toml`.

## 0.1.0 (2026-11-14)

Initial Python release (`ggtree-py`, imported as `ggtree`), ported from Bioconductor `ggtree` (`4.3.2`).

### Highlights
- **Tree parsers & structures (`ggtree.tree`)**: Newick (`read_tree`), BEAST NEXUS (`read_beast`), NHX (`read_nhx`), jplace (`read_jplace`), FASTA (`read_fasta`), `Phylo`, `TreeData`, `MultiPhylo`, `TreeDataList`, `Hclust`/`hclust`, and bundled sample datasets (`system_file`).
- **Tidy tree operations (`ggtree.tidytree`)**: `MRCA`, `groupOTU`, `groupClade`, `offspring`, `ancestor`, `parent`, `child`, `sibling`, `nodeid`, `nodelab`, `tree_subset`, `td_filter`, `td_mutate`, and `td_unnest`.
- **12+ tree layouts (`ggtree.layout`)**: `rectangular`, `slanted`, `circular`, `fan`, `inward_circular`, `dendrogram`, `equal_angle`, `daylight` (`unrooted`), `tidy`, `ellipse`, `roundrect`, `radial`, `ape`, and `tree_and_leaf`, plus time-scaled MRSD coordinates (`Date2decimal`, `decimal2Date`).
- **Grammar-of-graphics geoms & scales (`ggtree.aes`, `ggtree.geoms`, `ggtree.themes_scales`)**: Symbolic expressions (`aes`, `isTip`, `parent`), `geom_tree`, `geom_tiplab`, `geom_tiplab2`, `geom_nodelab`, `geom_point2`, `geom_tippoint`, `geom_nodepoint`, `geom_rootpoint`, `geom_rootedge`, `geom_aline`, `geom_hilight`, `geom_balance`, `geom_cladelab`, `geom_strip`, `geom_taxalink`, `geom_tanglelink`, `geom_range`, `geom_treescale`, `new_scale_fill`, `new_scale_color`, Viridis/Magma/Plasma/Inferno/Cividis/Brewer/gradient scales, and `theme_tree`, `theme_tree2`, `theme_dendrogram`.
- **Manipulation & composite plots (`ggtree.manipulate`, `ggtree.compose`)**: `collapse`, `expand`, `rotate`, `flip`, `scaleClade`, `viewClade`, `zoomClade`, `open_tree`, `rotate_tree`, `gheatmap`, `msaplot`, `geom_facet`/`facet_plot`, `nodepie`/`nodebar`/`inset`, `ggdensitree`, `as_tanglegram`/`ggdoubletree`, `multiplot`, `plot_list`, and `|`/`/` patchwork composition.
- **Zero-external-dependency SVG & PNG renderer (`ggtree.plot`)**: Direct vector SVG rendering and PNG export via `ggsave()`, `to_svg()`, and `to_png()` (using `rsvg-convert` when present, with pure-NumPy RGBA fallback).
- **Verified**: All 49 tests pass (`tests/test_examples.py` + `tests/test_ggtree.py`).
