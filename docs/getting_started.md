# Getting Started & Tree Data I/O

## 1. Reading & Creating Trees (`ggtree.tree`)

`ggtree` parses standard phylogenetic file formats and strings without external dependencies:

| Function | Format | Return Type |
| --- | --- | --- |
| `read_tree(file=None, text=None, ...)` | Newick / NEXUS trees | `Phylo` or `MultiPhylo` |
| `read_beast(file)` | BEAST annotated NEXUS (`[&height_95%_HPD=...]`, `posterior`, `rate`, ...) | `TreeData` or `TreeDataList` |
| `read_nhx(file=None, text=None)` | New Hampshire eXtended (`[&&NHX:...]`) | `TreeData` |
| `read_jplace(file)` | Evolutionary placement JSON (`.jplace`) | `TreeData` |
| `read_fasta(file)` | FASTA nucleotide / amino-acid alignment | `Dict[str, str]` |
| `write_tree(phy, file=None)` | Newick serializer | `str` |
| `rtree(n)` / `rmtree(N, n)` / `stree(n)` | Random & star/balanced trees (ported from `ape`) | `Phylo` / `MultiPhylo` |
| `hclust(d, method="complete", labels=None)` | Hierarchical clustering dendrograms (`as_phylo` supported) | `Hclust` |

Bundled sample datasets mirror `system.file(..., package="ggtree" / "treeio")`:

```python
import ggtree as gg

nwk = gg.system_file("extdata", "sample.nwk", package="treeio")
tr = gg.read_tree(nwk)

beast_tr = gg.read_beast(gg.system_file("extdata/MCC_FluA_H3.tree", package="ggtree"))
raxml_tr = gg.read_tree(gg.system_file("extdata/RAxML", "RAxML_bipartitions.H3", package="treeio"))
```

### Tree Inspection & Tidytree Conversion

- `Ntip(tr)`, `Nnode(tr)`, `Nnode2(tr)`, `is_rooted(tr)`, `get_tree(tr)`, `get_seqs(tr)`
- `as_phylo(x)`, `as_treedata(x)`, `as_tibble(x)` -> `TblTree` containing `parent`, `node`, `branch.length`, `label`, and any attached node annotations.
- `fortify(tr, layout="rectangular", ...)` -> `TblTree` with layout coordinates (`x`, `y`, `branch`, `angle`, `isTip`, ...).

---

## 2. Attaching Metadata & Transforming Tree Data

External annotations (such as sample metadata, traits, or host species) can be attached directly to a `TreeData` object or onto a `GGTree` plot:

### Operator `%` / `<<` / `@` (Port of R `%<+%`)

```python
import ggtree as gg

tr = gg.read_tree(text="((A:1,B:1):1,(C:1,D:1):1);")
meta = gg.DataFrame({
    "taxa": ["A", "B", "C", "D"],
    "score": [10, 20, 30, 40],
    "group": ["g1", "g1", "g2", "g2"],
})

p = (
    gg.ggtree(tr)
    % meta
    + gg.geom_tiplab(gg.aes(color="group"))
    + gg.geom_tippoint(gg.aes(size="score"))
)
```

- `p % meta` (or `p << meta`, `p @ meta`) joins `meta` onto `p.data` matching the label/taxa column (`taxa`, `label`, or first column).
- `p % new_tree` replaces the tree inside `p` (equivalent to R's `%<%`) while retaining layers and themes.

### Table Joins & Pipeable Verbs (`td_filter`, `td_mutate`, `td_unnest`)

- `left_join(tree, df, by="label")` / `full_join` / `inner_join` return an enriched `TreeData` (when joining onto `Phylo`/`TreeData`) or `DataFrame`.
- `td_filter(...)`, `td_mutate(...)`, and `td_unnest(cols)` create layer `data=` callables matching R's lazy `td_*()` helpers:

```python
p = gg.ggtree(tr) + gg.geom_point(data=gg.td_filter("isTip"), color="red")
```

---

## 3. Symbolic Aesthetics (`aes`) & Subsetting

`aes(...)` accepts:
- Symbolic column references imported from `ggtree` (`label`, `node`, `branch`, `x`, `y`, `angle`, `group`, `posterior`, `rate`, `height`, `value`, ...)
- Callable subset helpers (`isTip`, `parent`, negated with `~isTip` for `!isTip` internal nodes)
- String column names or expressions (e.g. `aes(color="group")`, `aes(subset="isTip & (node > 3)")`)

```python
p = (
    gg.ggtree(tr)
    + gg.geom_nodepoint(gg.aes(subset=~gg.isTip), color="orange", size=2)
    + gg.geom_text2(gg.aes(subset=~gg.isTip, label=gg.node), hjust=-0.3)
)
```

---

## 4. Building, Exporting SVG & PNG

- **SVG string**: `svg_text = p.to_svg(width=800, height=600)` (or `_repr_svg_()` for Jupyter notebooks)
- **PNG bytes**: `png_data = p.to_png(width=800, height=600)`
- **Save to disk (`ggsave`)**:
  ```python
  gg.ggsave("output.svg", plot=p, width=9.0, height=6.5, dpi=100)
  gg.ggsave("output.png", plot=p, width=9.0, height=6.5, dpi=100)
  gg.ggsave("output.html", plot=p, width=9.0, height=6.5)
  ```
  If no `plot` argument is supplied, `ggsave` uses `last_plot()`.
