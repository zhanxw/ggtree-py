"""
Phylogenetic tree data structures (Phylo, TreeData, MultiPhylo, TreeDataList, Hclust),
Newick/Nexus/BEAST/FASTA/jplace parsers, tree generators (rtree, rmtree, stree),
and tree manipulation utilities (ladderize, di2multi, as_polytomy, as_tibble, etc.).
"""

from __future__ import annotations

import copy
import json
import math
import os
import re
import warnings
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import numpy as np

from .dataframe import (
    ColumnSeries,
    DataFrame,
    Factor,
    TblTree,
    _is_scalar_na,
    _get_rng,
    c,
    full_join,
    left_join,
    tibble,
)


class EdgeArray(np.ndarray):
    """numpy.ndarray subclass for tr.edge that also exposes .length (for R's tr$edge.length)."""

    _owner: Optional["Phylo"] = None

    def __new__(cls, input_array: Any, owner: Optional["Phylo"] = None) -> "EdgeArray":
        obj = np.asarray(input_array, dtype=int).view(cls)
        obj._owner = owner
        return obj

    def __array_finalize__(self, obj: Any) -> None:
        if obj is None:
            return
        self._owner = getattr(obj, "_owner", None)

    @property
    def length(self) -> Optional[List[float]]:
        if self._owner is not None:
            return self._owner.edge_length
        return None

    @length.setter
    def length(self, val: Optional[Sequence[float]]) -> None:
        if self._owner is not None:
            self._owner.edge_length = val


class _LabelProxy:
    def __init__(self, owner: "Phylo", attr_name: str):
        self._owner = owner
        self._attr_name = attr_name

    @property
    def label(self) -> Any:
        return getattr(self._owner, self._attr_name)

    @label.setter
    def label(self, val: Any) -> None:
        setattr(self._owner, self._attr_name, val)


class _RootProxy:
    def __init__(self, owner: "Phylo"):
        self._owner = owner

    @property
    def edge(self) -> Optional[float]:
        return self._owner.root_edge

    @edge.setter
    def edge(self, val: Optional[float]) -> None:
        self._owner.root_edge = val


class Phylo:
    """
    R ape::phylo compatible phylogenetic tree object.
    1-based node numbering:
      - Tips: 1 .. Ntip
      - Root: Ntip + 1
      - Internal nodes: Ntip + 1 .. Ntip + Nnode
    Supports both Pythonic attributes (tip_label, edge_length, node_label, root_edge)
    and R-style dot/dict access (tr.tip.label, tr.edge.length, tr.root.edge, tr['tip.label']).
    """

    def __init__(
        self,
        edge: Any,
        tip_label: Sequence[str],
        Nnode: Optional[int] = None,
        edge_length: Optional[Sequence[Optional[float]]] = None,
        node_label: Optional[Sequence[Optional[str]]] = None,
        root_edge: Optional[float] = None,
        attrs: Optional[Dict[str, Any]] = None,
    ):
        edge_arr = np.asarray(edge, dtype=int)
        if edge_arr.ndim == 1:
            edge_arr = edge_arr.reshape(-1, 2)
        self._edge = EdgeArray(edge_arr, owner=self)
        self._tip_label = [str(x) for x in tip_label]
        if Nnode is None:
            internal_nodes = set(self._edge[:, 0].tolist()) | set(self._edge[:, 1].tolist())
            ntip = len(self._tip_label)
            self.Nnode = len([n for n in internal_nodes if n > ntip])
        else:
            self.Nnode = int(Nnode)
        self._edge_length: Optional[List[float]] = (
            [float(x) if x is not None and not _is_scalar_na(x) else float("nan") for x in edge_length]
            if edge_length is not None
            else None
        )
        self._node_label: Optional[List[Optional[str]]] = (
            [str(x) if x is not None and not _is_scalar_na(x) else None for x in node_label]
            if node_label is not None
            else None
        )
        self.root_edge: Optional[float] = float(root_edge) if root_edge is not None and not _is_scalar_na(root_edge) else None
        self.attrs: Dict[str, Any] = dict(attrs) if attrs else {}

    @property
    def edge(self) -> EdgeArray:
        self._edge._owner = self
        return self._edge

    @edge.setter
    def edge(self, val: Any) -> None:
        self._edge = EdgeArray(val, owner=self)

    @property
    def tip_label(self) -> List[str]:
        return self._tip_label

    @tip_label.setter
    def tip_label(self, val: Sequence[Any]) -> None:
        self._tip_label = [str(x) for x in val]

    @property
    def edge_length(self) -> Optional[List[float]]:
        return self._edge_length

    @edge_length.setter
    def edge_length(self, val: Optional[Sequence[Any]]) -> None:
        if val is None:
            self._edge_length = None
        else:
            self._edge_length = [float(x) if x is not None and not _is_scalar_na(x) else float("nan") for x in val]

    @property
    def node_label(self) -> Optional[List[Optional[str]]]:
        return self._node_label

    @node_label.setter
    def node_label(self, val: Optional[Sequence[Any]]) -> None:
        if val is None:
            self._node_label = None
        else:
            self._node_label = [str(x) if x is not None and not _is_scalar_na(x) else None for x in val]

    @property
    def tip(self) -> _LabelProxy:
        return _LabelProxy(self, "tip_label")

    @property
    def node(self) -> _LabelProxy:
        return _LabelProxy(self, "node_label")

    @property
    def root(self) -> _RootProxy:
        return _RootProxy(self)

    def copy(self) -> "Phylo":
        return Phylo(
            edge=self._edge.copy(),
            tip_label=list(self._tip_label),
            Nnode=self.Nnode,
            edge_length=list(self._edge_length) if self._edge_length is not None else None,
            node_label=list(self._node_label) if self._node_label is not None else None,
            root_edge=self.root_edge,
            attrs=copy.deepcopy(self.attrs),
        )

    def __getitem__(self, key: str) -> Any:
        k = key.replace(".", "_")
        if hasattr(self, k):
            return getattr(self, k)
        if key in self.attrs:
            return self.attrs[key]
        raise KeyError(key)

    def __setitem__(self, key: str, value: Any) -> None:
        k = key.replace(".", "_")
        if k in ("edge", "tip_label", "edge_length", "node_label", "root_edge", "Nnode"):
            setattr(self, k, value)
        else:
            self.attrs[key] = value

    def __repr__(self) -> str:
        ntip = len(self._tip_label)
        rooted_str = "rooted" if is_rooted(self) else "unrooted"
        br_str = "includes branch lengths" if self._edge_length is not None else "no branch lengths"
        return f"Phylogenetic tree with {ntip} tips and {self.Nnode} internal nodes ({rooted_str}; {br_str})."


class TreeData:
    """
    R treeio::treedata S4-compatible object combining a Phylo tree with node-level
    annotations (.data) and extra metadata (.extraInfo).
    """

    def __init__(
        self,
        phylo: Phylo,
        data: Optional[DataFrame] = None,
        extraInfo: Optional[DataFrame] = None,
        file: str = "",
        attrs: Optional[Dict[str, Any]] = None,
        tip_seq: Optional[Dict[str, str]] = None,
    ):
        self.phylo = phylo.copy() if isinstance(phylo, Phylo) else as_phylo(phylo)
        self.data = DataFrame(data) if data is not None else DataFrame()
        self.extraInfo = DataFrame(extraInfo) if extraInfo is not None else DataFrame()
        self.file = file
        self.attrs: Dict[str, Any] = dict(attrs) if attrs else {}
        self.tip_seq: Optional[Dict[str, str]] = dict(tip_seq) if tip_seq is not None else None

    @property
    def tip_label(self) -> List[str]:
        return self.phylo.tip_label

    @tip_label.setter
    def tip_label(self, val: Sequence[str]) -> None:
        self.phylo.tip_label = val

    @property
    def tip(self) -> _LabelProxy:
        return self.phylo.tip

    @property
    def node_label(self) -> Optional[List[Optional[str]]]:
        return self.phylo.node_label

    @node_label.setter
    def node_label(self, val: Optional[Sequence[Optional[str]]]) -> None:
        self.phylo.node_label = val

    @property
    def node(self) -> _LabelProxy:
        return self.phylo.node

    @property
    def root_edge(self) -> Optional[float]:
        return self.phylo.root_edge

    @root_edge.setter
    def root_edge(self, val: Optional[float]) -> None:
        self.phylo.root_edge = val

    @property
    def root(self) -> _RootProxy:
        return self.phylo.root

    @property
    def edge(self) -> EdgeArray:
        return self.phylo.edge

    @property
    def edge_length(self) -> Optional[List[float]]:
        return self.phylo.edge_length

    @property
    def Nnode(self) -> int:
        return self.phylo.Nnode

    def copy(self) -> "TreeData":
        return TreeData(
            phylo=self.phylo.copy(),
            data=self.data.copy(),
            extraInfo=self.extraInfo.copy(),
            file=self.file,
            attrs=copy.deepcopy(self.attrs),
            tip_seq=dict(self.tip_seq) if self.tip_seq is not None else None,
        )

    def __repr__(self) -> str:
        fields = [c for c in self.data.columns if c != "node"] + [c for c in self.extraInfo.columns if c != "node"]
        return f"'treedata' S4 object ({self.phylo!r}, associated data fields: {fields})."


class MultiPhylo(list):
    """R ape::multiPhylo list of Phylo trees."""

    def __init__(self, trees: Any = (), names: Optional[Sequence[str]] = None, **named_trees: Any):
        if isinstance(trees, dict):
            names = list(trees.keys()) if names is None else list(names)
            trees = list(trees.values())
        elif named_trees and not trees:
            names = list(named_trees.keys()) if names is None else list(names)
            trees = list(named_trees.values())
        super().__init__(trees)
        self.names: Optional[List[str]] = list(names) if names is not None else None
        self.attrs: Dict[str, Any] = {}


class TreeDataList(list):
    """R treeio::treedataList list of TreeData objects."""

    def __init__(self, trees: Any = (), names: Optional[Sequence[str]] = None, **named_trees: Any):
        if isinstance(trees, dict):
            names = list(trees.keys()) if names is None else list(names)
            trees = list(trees.values())
        elif named_trees and not trees:
            names = list(named_trees.keys()) if names is None else list(names)
            trees = list(named_trees.values())
        td_items = [as_treedata(t) if isinstance(t, Phylo) else t for t in trees]
        super().__init__(td_items)
        self.names: Optional[List[str]] = list(names) if names is not None else None
        self.attrs: Dict[str, Any] = {}


class Hclust:
    """
    R stats::hclust / dendrogram wrapper supporting conversion from
    scipy.cluster.hierarchy linkage matrices or distance matrices.
    """

    def __init__(
        self,
        merge: Any,
        height: Sequence[float],
        order: Optional[Sequence[int]] = None,
        labels: Optional[Sequence[str]] = None,
    ):
        self.merge = np.asarray(merge, dtype=int)
        self.height = [float(h) for h in height]
        n = len(self.height) + 1
        self.labels = [str(x) for x in labels] if labels is not None else [f"t{i+1}" for i in range(n)]
        self.order = list(order) if order is not None else list(range(1, n + 1))
        self.attrs: Dict[str, Any] = {}

    @classmethod
    def from_linkage(cls, Z: Any, labels: Optional[Sequence[str]] = None) -> "Hclust":
        """Create Hclust from a scipy-style (n-1, 4) linkage matrix Z."""
        Z_arr = np.asarray(Z, dtype=float)
        n = Z_arr.shape[0] + 1
        merge = np.zeros((n - 1, 2), dtype=int)
        height = []
        for i in range(n - 1):
            for j in (0, 1):
                idx = int(Z_arr[i, j])
                if idx < n:
                    merge[i, j] = -(idx + 1)
                else:
                    merge[i, j] = (idx - n) + 1
            height.append(float(Z_arr[i, 2]))
        return cls(merge=merge, height=height, labels=labels)

    def to_phylo(self, hang: float = -1) -> Phylo:
        """Convert Hclust to Phylo (matching ape::as.phylo.hclust + ggtree's fortify.hclust)."""
        N = len(self.height)
        ntip = N + 1
        edges: List[List[int]] = []
        edge_lengths: List[float] = []
        tip_edge_len = [0.0] * ntip

        # Node mapping: hclust step N is the root (node ntip + 1), step i is node ntip + N - i
        def step_to_node(step_1based: int) -> int:
            return ntip + (N - step_1based + 1)

        def build(step_1based: int) -> None:
            pnode = step_to_node(step_1based)
            p_height = self.height[step_1based - 1]
            for j in (0, 1):
                m = int(self.merge[step_1based - 1, j])
                if m < 0:
                    tip_id = -m
                    edges.append([pnode, tip_id])
                    # ape::as.phylo.hclust divides height by 2 (or uses height difference)
                    edge_lengths.append(p_height / 2.0)
                    if hang > 0:
                        tip_edge_len[tip_id - 1] = max(0.0, (p_height - hang * max(self.height)) / 2.0)
                else:
                    cnode = step_to_node(m)
                    c_height = self.height[m - 1]
                    edges.append([pnode, cnode])
                    edge_lengths.append((p_height - c_height) / 2.0)
                    build(m)

        build(N)
        phy = Phylo(edge=edges, tip_label=self.labels, Nnode=N, edge_length=edge_lengths)
        if hang > 0:
            phy.attrs["tip.edge.len"] = tip_edge_len
        return phy


def hclust(d: Any, method: str = "complete", labels: Optional[Sequence[str]] = None) -> Hclust:
    """Simple agglomerative hierarchical clustering (like R's stats::hclust)."""
    mat = np.asarray(d, dtype=float)
    if mat.ndim == 1:
        # Condensed distance vector of length n*(n-1)//2
        m = len(mat)
        n = int((1 + math.isqrt(1 + 8 * m)) // 2)
        full = np.zeros((n, n), dtype=float)
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                full[i, j] = full[j, i] = mat[k]
                k += 1
        mat = full
    n = mat.shape[0]
    clusters: Dict[int, List[int]] = {-i: [i - 1] for i in range(1, n + 1)}
    merge = np.zeros((n - 1, 2), dtype=int)
    heights: List[float] = []

    for step in range(1, n):
        keys = list(clusters.keys())
        best_d = float("inf")
        best_pair = (keys[0], keys[1])
        for idx1 in range(len(keys)):
            for idx2 in range(idx1 + 1, len(keys)):
                k1, k2 = keys[idx1], keys[idx2]
                pts1, pts2 = clusters[k1], clusters[k2]
                sub = mat[np.ix_(pts1, pts2)]
                if method == "single":
                    dist = float(np.min(sub))
                elif method in ("average", "upgma"):
                    dist = float(np.mean(sub))
                else:
                    dist = float(np.max(sub))
                if dist < best_d:
                    best_d = dist
                    best_pair = (k1, k2)
        k1, k2 = best_pair
        merge[step - 1, 0] = k1
        merge[step - 1, 1] = k2
        heights.append(best_d)
        clusters[step] = clusters.pop(k1) + clusters.pop(k2)

    return Hclust(merge=merge, height=heights, labels=labels)


# Tree inspection & conversion functions
def Ntip(x: Any) -> int:
    """Return number of tips in a tree object."""
    if isinstance(x, Phylo):
        return len(x.tip_label)
    if isinstance(x, TreeData):
        return len(x.phylo.tip_label)
    if isinstance(x, (MultiPhylo, TreeDataList, list)):
        return [Ntip(item) for item in x]  # type: ignore[return-value]
    if isinstance(x, DataFrame):
        if "isTip" in x:
            return int(x["isTip"].sum())
        parents = set(x["parent"].tolist())
        return sum(1 for n in x["node"].tolist() if n not in parents)
    if isinstance(x, Hclust):
        return len(x.labels)
    if hasattr(x, "data") and isinstance(x.data, DataFrame):
        return Ntip(x.data)
    raise TypeError(f"Cannot determine Ntip for {type(x)}")


def Nnode(x: Any, internal_only: bool = True, **kwargs: Any) -> int:
    """Return number of internal nodes (or total nodes if internal.only=False)."""
    internal_only = kwargs.get("internal.only", internal_only)
    if isinstance(x, Phylo):
        return x.Nnode if internal_only else len(x.tip_label) + x.Nnode
    if isinstance(x, TreeData):
        return Nnode(x.phylo, internal_only=internal_only)
    if isinstance(x, DataFrame):
        total = len(x)
        if not internal_only:
            return total
        return total - Ntip(x)
    if isinstance(x, Hclust):
        return len(x.height) if internal_only else len(x.labels) + len(x.height)
    if hasattr(x, "data") and isinstance(x.data, DataFrame):
        return Nnode(x.data, internal_only=internal_only)
    raise TypeError(f"Cannot determine Nnode for {type(x)}")


def Nnode2(x: Any) -> int:
    """Return total number of nodes (tips + internal nodes), matching treeio::Nnode2."""
    return Nnode(x, internal_only=False)


def getNodeNum(x: Any) -> int:
    """Return total number of nodes (tips + internal nodes)."""
    return Nnode(x, internal_only=False)


def is_rooted(x: Any) -> bool:
    """Check if tree is rooted."""
    phy = as_phylo(x)
    root = len(phy.tip_label) + 1
    root_children = np.sum(phy.edge[:, 0] == root)
    return bool(root_children <= 2 or phy.root_edge is not None)


def get_tree(x: Any) -> Phylo:
    """Extract Phylo object from Phylo, TreeData, or DataFrame."""
    return as_phylo(x)


def as_phylo(x: Any, **kwargs: Any) -> Phylo:
    """Convert Phylo, TreeData, DataFrame/TblTree, Hclust, or Newick string to Phylo."""
    if isinstance(x, Phylo):
        return x
    if isinstance(x, TreeData):
        return x.phylo
    if isinstance(x, Hclust):
        return x.to_phylo(hang=kwargs.get("hang", -1))
    if isinstance(x, str):
        return read_tree(text=x)
    if hasattr(x, "data") and isinstance(x.data, DataFrame):
        x = x.data
    if isinstance(x, DataFrame):
        # Convert tbl_tree DataFrame back to Phylo (matching tidytree::as.phylo.tbl_tree)
        df = x.sort_values("node")
        nodes = [int(v) for v in df["node"].tolist()]
        parents = [int(v) for v in df["parent"].tolist()]
        labels = df["label"].tolist() if "label" in df else [f"t{i}" for i in nodes]
        if "isTip" in df:
            is_tip = [bool(v) for v in df["isTip"].tolist()]
        else:
            pset = set(parents)
            is_tip = [n not in pset or (p == n and False) for p, n in zip(parents, nodes)]
        tip_nodes = [n for n, flag in zip(nodes, is_tip) if flag]
        ntip = len(tip_nodes)
        tip_labels = [str(labels[nodes.index(n)]) for n in sorted(tip_nodes)]
        internal_nodes = sorted([n for n, flag in zip(nodes, is_tip) if not flag])
        node_labels = [
            str(labels[nodes.index(n)]) if not _is_scalar_na(labels[nodes.index(n)]) else None
            for n in internal_nodes
        ]
        edges = []
        edge_lens = []
        has_brlen = "branch.length" in df
        br_col = df["branch.length"].tolist() if has_brlen else []
        root_edge = None
        for idx, (p, n) in enumerate(zip(parents, nodes)):
            if p == n:
                if has_brlen and not _is_scalar_na(br_col[idx]):
                    root_edge = float(br_col[idx])
                continue
            edges.append([p, n])
            if has_brlen:
                val = br_col[idx]
                edge_lens.append(float(val) if not _is_scalar_na(val) else 0.0)
        return Phylo(
            edge=edges,
            tip_label=tip_labels,
            Nnode=len(internal_nodes),
            edge_length=edge_lens if has_brlen else None,
            node_label=node_labels if any(x is not None for x in node_labels) else None,
            root_edge=root_edge,
        )
    # Support networkx / igraph graphs if passed
    if hasattr(x, "edges") and callable(getattr(x, "edges", None)):
        edges_raw = list(x.edges())
        all_v = []
        for u, v in edges_raw:
            if u not in all_v:
                all_v.append(u)
            if v not in all_v:
                all_v.append(v)
        parents_v = {u for u, _ in edges_raw}
        tips_v = [v for v in all_v if v not in parents_v]
        internals_v = [v for v in all_v if v in parents_v]
        v_map = {}
        for i, tv in enumerate(tips_v, 1):
            v_map[tv] = i
        for i, iv in enumerate(internals_v, len(tips_v) + 1):
            v_map[iv] = i
        edge_mat = [[v_map[u], v_map[v]] for u, v in edges_raw]
        return Phylo(edge=edge_mat, tip_label=[str(tv) for tv in tips_v], Nnode=len(internals_v))
    raise TypeError(f"Cannot convert {type(x)} to Phylo")


def as_treedata(x: Any, **kwargs: Any) -> TreeData:
    """Convert Phylo, TreeData, or TblTree to TreeData."""
    if isinstance(x, TreeData):
        return x.copy()
    if isinstance(x, Phylo):
        return TreeData(phylo=x.copy())
    if isinstance(x, DataFrame):
        phy = as_phylo(x)
        core_cols = {"parent", "node", "branch.length", "label", "isTip", "x", "y", "branch", "angle"}
        extra_cols = ["node"] + [c for c in x.columns if c not in core_cols]
        data = x._select_cols(extra_cols) if len(extra_cols) > 1 else DataFrame()
        return TreeData(phylo=phy, data=data)
    raise TypeError(f"Cannot convert {type(x)} to TreeData")


def as_tibble(x: Any, **kwargs: Any) -> TblTree:
    """
    Convert Phylo or TreeData to TblTree (exact port of tidytree::as_tibble.phylo
    and tidytree::as_tibble.treedata).
    """
    if isinstance(x, TblTree):
        return x.copy()
    if isinstance(x, DataFrame):
        return TblTree(x, attrs=x.attrs)
    if isinstance(x, TreeData):
        res = as_tibble(x.phylo)
        for df_anno in (x.data, x.extraInfo):
            if df_anno is not None and len(df_anno) > 0 and "node" in df_anno:
                by = ["node", "parent"] if "parent" in df_anno.columns else "node"
                res = full_join(res, df_anno, by=by)
        return TblTree(res, attrs=res.attrs)
    if isinstance(x, Phylo):
        phylo = x
        ntip = Ntip(phylo)
        N = Nnode(phylo, internal_only=False)
        edge = phylo.edge
        parent_map: Dict[int, int] = {}
        brlen_map: Dict[int, Optional[float]] = {}
        has_brlen = phylo.edge_length is not None
        for idx in range(edge.shape[0]):
            p = int(edge[idx, 0])
            ch = int(edge[idx, 1])
            parent_map[ch] = p
            if has_brlen and phylo.edge_length is not None:
                brlen_map[ch] = float(phylo.edge_length[idx])

        labels: List[Optional[str]] = [None] * N
        for i in range(min(ntip, len(phylo.tip_label))):
            labels[i] = phylo.tip_label[i]
        if phylo.node_label is not None:
            for j, nl in enumerate(phylo.node_label):
                if ntip + j < N:
                    labels[ntip + j] = nl

        nodes_col = list(range(1, N + 1))
        parents_col = [parent_map.get(nd, nd) for nd in nodes_col]
        cols_dict: Dict[str, Any] = {
            "parent": parents_col,
            "node": nodes_col,
        }
        if has_brlen:
            br_col: List[Optional[float]] = []
            for nd, p in zip(nodes_col, parents_col):
                if nd == p:
                    br_col.append(phylo.root_edge if phylo.root_edge is not None else None)
                else:
                    br_col.append(brlen_map.get(nd, None))
            cols_dict["branch.length"] = br_col
        cols_dict["label"] = labels

        res = TblTree(cols_dict)
        # Attach custom node-level attributes stored on phylo.attrs (like groupOTU / groupClade)
        reserved = {"names", "class", "order", "reroot", "node_map", "tip.edge.len"}
        for attr_name, attr_val in phylo.attrs.items():
            if attr_name not in reserved and isinstance(attr_val, (list, tuple, np.ndarray, Factor, ColumnSeries)):
                if len(attr_val) == N:
                    res[attr_name] = attr_val
        return res
    raise TypeError(f"Cannot convert {type(x)} to TblTree")


# Reorder & Ladderize
def reorder_phylo(phy: Phylo, order: str = "cladewise") -> Phylo:
    """
    Reorder rows of phy.edge in 'cladewise' (pre-order) or 'postorder'
    (matching ape::reorder.phylo).
    """
    out = phy.copy()
    edge = out.edge
    has_len = out.edge_length is not None
    len_map: Dict[Tuple[int, int], float] = {}
    if has_len and out.edge_length is not None:
        for i in range(edge.shape[0]):
            len_map[(int(edge[i, 0]), int(edge[i, 1]))] = float(out.edge_length[i])

    children_of: Dict[int, List[int]] = {}
    all_parents = set()
    all_children = set()
    for i in range(edge.shape[0]):
        p, c_node = int(edge[i, 0]), int(edge[i, 1])
        children_of.setdefault(p, []).append(c_node)
        all_parents.add(p)
        all_children.add(c_node)

    roots = list(all_parents - all_children)
    root = roots[0] if roots else len(out.tip_label) + 1

    new_edges: List[List[int]] = []
    if order == "cladewise":
        # Iterative pre-order DFS
        stack = [root]
        visited = set()
        while stack:
            u = stack.pop()
            if u in visited:
                continue
            visited.add(u)
            chs = children_of.get(u, [])
            # Push in reverse so first child is popped first, but record edges as traversed
            for ch in reversed(chs):
                stack.append((u, ch))  # type: ignore[arg-type]
        # Let's do clean edge stack for cladewise
        edge_stack: List[Tuple[int, int]] = [(root, ch) for ch in reversed(children_of.get(root, []))]
        while edge_stack:
            p, c_node = edge_stack.pop()
            new_edges.append([p, c_node])
            for gc in reversed(children_of.get(c_node, [])):
                edge_stack.append((c_node, gc))
    elif order == "postorder":
        # Iterative post-order traversal of internal nodes, emitting children edges of each node
        stack_nodes = [(root, False)]
        while stack_nodes:
            u, done = stack_nodes.pop()
            if done:
                for ch in children_of.get(u, []):
                    new_edges.append([u, ch])
            else:
                stack_nodes.append((u, True))
                for ch in reversed(children_of.get(u, [])):
                    if ch in children_of:
                        stack_nodes.append((ch, False))
    else:
        raise ValueError(f"Unknown order: {order}")

    if len(new_edges) == edge.shape[0]:
        out.edge = new_edges
        if has_len:
            out.edge_length = [len_map[(p, c_node)] for p, c_node in new_edges]
    out.attrs["order"] = order
    return out


def ladderize(phy: Phylo, right: bool = True) -> Phylo:
    """
    Exact port of ape::ladderize(phy, right = TRUE/FALSE).
    """
    out = reorder_phylo(phy, "postorder")
    edge = out.edge
    ntip = len(out.tip_label)
    nnode_total = int(np.max(edge)) if edge.size > 0 else ntip

    # Compute N[node] = number of descendant tips (matching ape's node_depth method=1)
    N = np.zeros(nnode_total + 1, dtype=float)
    N[1 : ntip + 1] = 1.0
    for i in range(edge.shape[0]):
        p, ch = int(edge[i, 0]), int(edge[i, 1])
        N[p] += N[ch]

    # Sort edges by (parent, N[child]) with decreasing = right (matching R's order(x, y, decreasing=right))
    indices = list(range(edge.shape[0]))
    if right:
        indices.sort(key=lambda idx: (-int(edge[idx, 0]), -N[int(edge[idx, 1])]))
    else:
        indices.sort(key=lambda idx: (int(edge[idx, 0]), N[int(edge[idx, 1])]))

    desc: Dict[int, List[int]] = {}
    for idx in indices:
        p, ch = int(edge[idx, 0]), int(edge[idx, 1])
        desc.setdefault(p, []).append(ch)

    el_map: Dict[int, float] = {}
    if out.edge_length is not None:
        for i in range(edge.shape[0]):
            el_map[int(edge[i, 1])] = float(out.edge_length[i])

    # Cladewise reconstruction from root (ntip + 1)
    root = ntip + 1
    if root not in desc and desc:
        all_p = set(desc.keys())
        all_c = {c_node for chs in desc.values() for c_node in chs}
        roots = list(all_p - all_c)
        if roots:
            root = roots[0]

    new_edges: List[List[int]] = []
    edge_stack: List[Tuple[int, int]] = [(root, ch) for ch in reversed(desc.get(root, []))]
    while edge_stack:
        p, ch = edge_stack.pop()
        new_edges.append([p, ch])
        for gc in reversed(desc.get(ch, [])):
            edge_stack.append((ch, gc))

    if len(new_edges) == edge.shape[0]:
        out.edge = new_edges
        if out.edge_length is not None:
            out.edge_length = [el_map[ch] for _, ch in new_edges]
    out.attrs["order"] = "cladewise"
    return out


# Newick / Nexus / BEAST Parser
def _strip_comments_preserve_beast(text: str) -> str:
    """Keep [&...] extended Newick comments while stripping ordinary [...] comments."""
    return re.sub(r"\[(?!!|&)[^\]]*\]", "", text)


def _parse_beast_annotation(comment: str) -> Dict[str, Any]:
    """Parse [&key=val,key2={v1,v2},...] BEAST/NHX comment into a dict."""
    inner = comment.strip()
    if inner.startswith("["):
        inner = inner[1:]
    if inner.endswith("]"):
        inner = inner[:-1]
    if inner.startswith("&") or inner.startswith("!"):
        inner = inner[1:]
    if inner.startswith("NHX:"):
        inner = inner[4:]
        sep_char = ":"
    else:
        sep_char = ","

    result: Dict[str, Any] = {}
    # Split by top-level commas (outside {...} and "...")
    parts: List[str] = []
    cur: List[str] = []
    brace_depth = 0
    in_quote = False
    for ch in inner:
        if ch == '"':
            in_quote = not in_quote
            cur.append(ch)
        elif not in_quote and ch == "{":
            brace_depth += 1
            cur.append(ch)
        elif not in_quote and ch == "}":
            brace_depth = max(0, brace_depth - 1)
            cur.append(ch)
        elif not in_quote and brace_depth == 0 and ch == sep_char:
            parts.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    if cur:
        parts.append("".join(cur).strip())

    for part in parts:
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        k = k.strip()
        v = v.strip()
        if v.startswith("{") and v.endswith("}"):
            items = [x.strip().strip('"').strip("'") for x in v[1:-1].split(",") if x.strip()]
            parsed_items: List[Any] = []
            for it in items:
                try:
                    parsed_items.append(float(it))
                except ValueError:
                    parsed_items.append(it)
            result[k] = parsed_items
        else:
            v_clean = v.strip('"').strip("'")
            try:
                if "." in v_clean or "e" in v_clean.lower():
                    result[k] = float(v_clean)
                else:
                    result[k] = int(v_clean)
            except ValueError:
                result[k] = v_clean
    return result


def _parse_single_newick(
    nwk: str,
    translate: Optional[Dict[str, str]] = None,
    extract_annotations: bool = False,
) -> Union[Phylo, TreeData]:
    """
    Iterative Newick / Extended Newick parser that matches ape::read.tree's exact
    1-based node numbering (tips 1..Ntip in left-to-right order, internal nodes
    Ntip+1..Ntip+Nnode in pre-order as '(' is opened, edges in pre-order).
    """
    s = nwk.strip()
    if s.endswith(";"):
        s = s[:-1].strip()

    # Tokenize Newick string with [&...] comments
    token_re = re.compile(
        r"""
        (\[&[^\]]*\]|\[![^\]]*\]|\[[^\]]*\]) |  # comment / annotation
        ([(),:;])                             |  # punctuation
        ('[^']*'|"[^"]*")                     |  # quoted label
        ([^\s(),:;\[\]]+)                        # unquoted token
        """,
        re.VERBOSE,
    )
    raw_tokens = [m.group(0) for m in token_re.finditer(s)]

    # Build AST iteratively
    # Each AST node: {"children": [], "label": None, "length": None, "ann": {}, "pre_idx": int}
    root_ast: Dict[str, Any] = {"children": [], "label": None, "length": None, "ann": {}, "pre_idx": 0}
    stack: List[Dict[str, Any]] = [root_ast]
    current = root_ast
    internal_counter = 0
    tips_list: List[Dict[str, Any]] = []
    internals_list: List[Dict[str, Any]] = []
    has_any_length = False

    i = 0
    n_tok = len(raw_tokens)
    while i < n_tok:
        tok = raw_tokens[i]
        if tok == "(":
            if current is root_ast and internal_counter == 0 and len(root_ast["children"]) == 0:
                current["pre_idx"] = internal_counter
                internal_counter += 1
                internals_list.append(current)
            else:
                new_node: Dict[str, Any] = {
                    "children": [],
                    "label": None,
                    "length": None,
                    "ann": {},
                    "pre_idx": internal_counter,
                }
                internal_counter += 1
                internals_list.append(new_node)
                current["children"].append(new_node)
                stack.append(new_node)
                current = new_node
            i += 1
        elif tok == ",":
            if len(stack) > 1 and current is not stack[-1]:
                pass
            current = stack[-1]
            # Wait: after a leaf child or closed internal child, current is the parent on top of stack
            i += 1
        elif tok == ")":
            closed = stack.pop() if len(stack) > 1 else root_ast
            current = closed
            # Check for optional post-parenthesis label, comment, or :branch_length
            i += 1
            while i < n_tok and raw_tokens[i] not in (",", ")", "("):
                sub = raw_tokens[i]
                if sub.startswith("["):
                    current["ann"].update(_parse_beast_annotation(sub))
                    i += 1
                elif sub == ":":
                    i += 1
                    while i < n_tok and raw_tokens[i].startswith("["):
                        current["ann"].update(_parse_beast_annotation(raw_tokens[i]))
                        i += 1
                    if i < n_tok:
                        try:
                            current["length"] = float(raw_tokens[i])
                            has_any_length = True
                        except ValueError:
                            pass
                        i += 1
                else:
                    lbl = sub.strip("'").strip('"')
                    current["label"] = translate.get(lbl, lbl) if translate else lbl
                    i += 1
            if stack:
                current = stack[-1]
        elif tok.startswith("["):
            current["ann"].update(_parse_beast_annotation(tok))
            i += 1
        else:
            # Leaf token!
            lbl = tok.strip("'").strip('"')
            if translate and lbl in translate:
                lbl = translate[lbl]
            leaf: Dict[str, Any] = {
                "children": [],
                "label": lbl,
                "length": None,
                "ann": {},
                "tip_idx": len(tips_list) + 1,
            }
            tips_list.append(leaf)
            current["children"].append(leaf)
            i += 1
            while i < n_tok and raw_tokens[i] not in (",", ")", "("):
                sub = raw_tokens[i]
                if sub.startswith("["):
                    leaf["ann"].update(_parse_beast_annotation(sub))
                    i += 1
                elif sub == ":":
                    i += 1
                    while i < n_tok and raw_tokens[i].startswith("["):
                        leaf["ann"].update(_parse_beast_annotation(raw_tokens[i]))
                        i += 1
                    if i < n_tok:
                        try:
                            leaf["length"] = float(raw_tokens[i])
                            has_any_length = True
                        except ValueError:
                            pass
                        i += 1
                else:
                    i += 1

    ntip = len(tips_list)
    nnode = len(internals_list)
    for int_node in internals_list:
        int_node["node_id"] = ntip + 1 + int_node["pre_idx"]
    for tip_node in tips_list:
        tip_node["node_id"] = tip_node["tip_idx"]

    # Traverse in pre-order to build edge matrix and edge_length
    edges: List[List[int]] = []
    edge_lengths: List[float] = []
    ann_records: List[Dict[str, Any]] = []

    work_stack: List[Dict[str, Any]] = [root_ast]
    while work_stack:
        u = work_stack.pop()
        u_id = u["node_id"]
        if u["ann"]:
            rec = dict(u["ann"])
            rec["node"] = u_id
            ann_records.append(rec)
        chs = u["children"]
        for ch in chs:
            edges.append([u_id, ch["node_id"]])
            edge_lengths.append(float(ch["length"]) if ch["length"] is not None else 0.0)
        for ch in reversed(chs):
            if ch["children"]:
                work_stack.append(ch)
            elif ch["ann"]:
                rec = dict(ch["ann"])
                rec["node"] = ch["node_id"]
                ann_records.append(rec)

    tip_labels = [str(t["label"]) for t in tips_list]
    node_labels = [int_n["label"] for int_n in internals_list]
    has_node_labels = any(nl is not None for nl in node_labels)
    root_edge = float(root_ast["length"]) if root_ast["length"] is not None else None

    phy = Phylo(
        edge=edges,
        tip_label=tip_labels,
        Nnode=nnode,
        edge_length=edge_lengths if has_any_length else None,
        node_label=node_labels if has_node_labels else None,
        root_edge=root_edge,
    )

    if extract_annotations:
        all_ann_keys: List[str] = ["node"]
        for r in ann_records:
            for k in r.keys():
                if k not in all_ann_keys:
                    all_ann_keys.append(k)
        df_ann = DataFrame(ann_records, columns=all_ann_keys) if ann_records else DataFrame({"node": list(range(1, ntip + nnode + 1))})
        return TreeData(phylo=phy, data=df_ann)
    return phy


def read_tree(
    file: Optional[str] = None,
    text: Optional[str] = None,
    tree_names: Optional[Sequence[str]] = None,
    keep_multi: bool = False,
    **kwargs: Any,
) -> Union[Phylo, MultiPhylo]:
    """
    Read Newick or Nexus phylogenetic tree(s) from file or text string
    (compatible with ape::read.tree and treeio::read.tree).
    """
    if text is None:
        if file is None:
            raise ValueError("Either 'file' or 'text' must be provided to read_tree()")
        with open(file, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
    else:
        content = text

    content_stripped = content.strip()
    if content_stripped.upper().startswith("#NEXUS"):
        trees = _parse_nexus_content(content_stripped, extract_annotations=False)
        if len(trees) == 1 and not keep_multi:
            return trees[0]  # type: ignore[return-value]
        return MultiPhylo(trees)

    # Strip non-annotation comments
    cleaned = re.sub(r"\[[^\]]*\]", "", content_stripped)
    # Split multiple Newick trees by ';'
    parts = [p.strip() for p in cleaned.split(";") if p.strip()]
    parsed = [_parse_single_newick(p + ";", extract_annotations=False) for p in parts]
    if len(parsed) == 1 and not keep_multi:
        return parsed[0]  # type: ignore[return-value]
    return MultiPhylo(parsed, names=tree_names)


def _parse_nexus_content(content: str, extract_annotations: bool = True) -> List[Union[Phylo, TreeData]]:
    translate: Dict[str, str] = {}
    trans_match = re.search(r"TRANSLATE\s+(.*?);", content, re.IGNORECASE | re.DOTALL)
    if trans_match:
        trans_body = trans_match.group(1)
        # Remove comments inside TRANSLATE
        trans_body = re.sub(r"\[[^\]]*\]", "", trans_body)
        # Split entries on comma or newline outside single/double quotes (NEXUS
        # Translate blocks may be comma-separated or newline-separated, e.g. BEAST)
        segments: List[str] = []
        buf: List[str] = []
        in_quote: Optional[str] = None
        for ch in trans_body:
            if ch in ("'", '"'):
                if in_quote is None:
                    in_quote = ch
                elif in_quote == ch:
                    in_quote = None
                buf.append(ch)
            elif in_quote is None and ch in (",", "\n", "\r"):
                seg = "".join(buf).strip()
                if seg:
                    segments.append(seg)
                buf.clear()
            else:
                buf.append(ch)
        if buf:
            seg = "".join(buf).strip()
            if seg:
                segments.append(seg)
        for item in segments:
            tokens = item.split(None, 1)
            if len(tokens) == 2:
                k = tokens[0].strip().rstrip(",").strip("'").strip('"')
                v = tokens[1].strip().rstrip(",").strip("'").strip('"')
                translate[k] = v

    tree_matches = re.finditer(
        r"TREE\s+\*?\s*([^\s=]+)\s*=\s*(?:\[[^\]]*\]\s*)*(\(.*?;)",
        content,
        re.IGNORECASE | re.DOTALL,
    )
    results: List[Union[Phylo, TreeData]] = []
    for m in tree_matches:
        nwk_str = m.group(2)
        results.append(_parse_single_newick(nwk_str, translate=translate, extract_annotations=extract_annotations))
    return results


def read_beast(file: str) -> Union[TreeData, TreeDataList]:
    """Read BEAST Nexus tree file into a TreeData object (compatible with treeio::read.beast)."""
    with open(file, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()
    trees = _parse_nexus_content(content, extract_annotations=True)
    for t in trees:
        if isinstance(t, TreeData):
            t.file = file
    if len(trees) == 1:
        return trees[0]  # type: ignore[return-value]
    return TreeDataList(trees)


def read_nhx(file: Optional[str] = None, text: Optional[str] = None) -> TreeData:
    """Read New Hampshire eXtended (NHX) tree into a TreeData object."""
    if text is None:
        if file is None:
            raise ValueError("Provide file or text")
        with open(file, "r", encoding="utf-8", errors="replace") as f:
            text = f.read()
    td = _parse_single_newick(text, extract_annotations=True)
    assert isinstance(td, TreeData)
    if file:
        td.file = file
    return td


def read_jplace(file: str) -> TreeData:
    """Read a .jplace phylogenetic placement file into a TreeData object."""
    with open(file, "r", encoding="utf-8") as f:
        obj = json.load(f)
    tree_str = obj.get("tree", "")
    # Replace {edge_num} with [&&NHX:edge_num=...]
    cleaned = re.sub(r"\{(\d+)\}", r"[&edge_num=\1]", tree_str)
    td = _parse_single_newick(cleaned, extract_annotations=True)
    assert isinstance(td, TreeData)
    td.file = file
    return td


def read_fasta(file: str) -> Dict[str, str]:
    """Read a FASTA alignment file into a dict mapping sequence ID -> sequence string."""
    seqs: Dict[str, str] = {}
    cur_id: Optional[str] = None
    cur_chunks: List[str] = []
    with open(file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_chunks)
                cur_id = line[1:].split()[0]
                cur_chunks = []
            else:
                cur_chunks.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(cur_chunks)
    return seqs


def write_tree(phy: Any, file: Optional[str] = None) -> str:
    """Export a Phylo or TreeData object to a Newick string."""
    p = as_phylo(phy)
    ntip = len(p.tip_label)
    children_of: Dict[int, List[int]] = {}
    br_of: Dict[int, Optional[float]] = {}
    edge = p.edge
    for i in range(edge.shape[0]):
        u, v = int(edge[i, 0]), int(edge[i, 1])
        children_of.setdefault(u, []).append(v)
        if p.edge_length is not None:
            br_of[v] = float(p.edge_length[i])

    root = ntip + 1

    def fmt_node(nd: int) -> str:
        if nd <= ntip:
            name = p.tip_label[nd - 1]
        else:
            chs = children_of.get(nd, [])
            inner = ",".join(fmt_node(c_node) for c_node in chs)
            lbl = ""
            if p.node_label is not None and (nd - ntip - 1) < len(p.node_label):
                nl = p.node_label[nd - ntip - 1]
                if nl is not None:
                    lbl = str(nl)
            name = f"({inner}){lbl}"
        if nd in br_of and br_of[nd] is not None:
            return f"{name}:{br_of[nd]:g}"
        return name

    nwk = fmt_node(root)
    if p.root_edge is not None:
        nwk += f":{p.root_edge:g}"
    nwk += ";"
    if file is not None:
        with open(file, "w", encoding="utf-8") as f:
            f.write(nwk + "\n")
    return nwk


# Random & Structured Tree Generators (ape::rtree, ape::rmtree, ape::stree)
def rtree(
    n: int,
    rooted: bool = True,
    tip_label: Optional[Sequence[str]] = None,
    br: Union[Callable[[int], Sequence[float]], str, None] = "runif",
    **kwargs: Any,
) -> Phylo:
    """
    Generate a random phylogenetic tree with n tips (compatible with ape::rtree).
    """
    tip_label = kwargs.get("tip.label", tip_label)
    rng = _get_rng()
    n = int(n)
    if n < 2:
        raise ValueError("Tree must have at least 2 tips")
    if tip_label is None:
        labels = [f"t{i}" for i in range(1, n + 1)]
        # Shuffle tip labels randomly like ape::rtree
        perm = rng.permutation(n)
        labels = [labels[int(i)] for i in perm]
    else:
        labels = [str(x) for x in tip_label]

    # Build random binary topology by recursively splitting n tips into (a, n-a)
    edges: List[List[int]] = []
    next_tip = 1
    next_internal = n + 1

    def split_clade(k: int) -> int:
        nonlocal next_tip, next_internal
        if k == 1:
            t_id = next_tip
            next_tip += 1
            return t_id
        my_id = next_internal
        next_internal += 1
        left_k = int(rng.integers(1, k))
        right_k = k - left_k
        c1 = split_clade(left_k)
        c2 = split_clade(right_k)
        edges.append([my_id, c1])
        edges.append([my_id, c2])
        return my_id

    split_clade(n)
    n_edges = len(edges)
    if br is None:
        edge_len = None
    elif callable(br):
        edge_len = [float(x) for x in br(n_edges)]
    else:
        edge_len = rng.uniform(0.0, 1.0, size=n_edges).tolist()

    phy = Phylo(edge=edges, tip_label=labels, Nnode=n - 1, edge_length=edge_len)
    return reorder_phylo(phy, "cladewise")


def rmtree(N: int, n: int, rooted: bool = True, tip_label: Optional[Sequence[str]] = None, **kwargs: Any) -> MultiPhylo:
    """Generate N random trees with n tips each (compatible with ape::rmtree)."""
    return MultiPhylo([rtree(n, rooted=rooted, tip_label=tip_label, **kwargs) for _ in range(int(N))])


def stree(n: int, type: str = "left", tip_label: Optional[Sequence[str]] = None, **kwargs: Any) -> Phylo:
    """
    Generate a regular tree ('left', 'right', 'balanced', 'star') with n tips
    (compatible with ape::stree).Iterative implementation supports n >= 2000+.
    """
    tip_label = kwargs.get("tip.label", tip_label)
    n = int(n)
    labels = [str(x) for x in tip_label] if tip_label is not None else [f"t{i}" for i in range(1, n + 1)]

    if type == "star":
        root = n + 1
        edges = [[root, i] for i in range(1, n + 1)]
        return Phylo(edge=edges, tip_label=labels, Nnode=1)
    elif type == "left":
        # Iterative caterpillar tree: root (n+1) -> tip n and internal (n+2), etc.
        edges = []
        for k in range(1, n - 1):
            pnode = n + k
            edges.append([pnode, pnode + 1])
            edges.append([pnode, n - k + 1])
        edges.append([2 * n - 1, 1])
        edges.append([2 * n - 1, 2])
        return Phylo(edge=edges, tip_label=labels, Nnode=n - 1)
    elif type == "right":
        edges = []
        for k in range(1, n - 1):
            pnode = n + k
            edges.append([pnode, k])
            edges.append([pnode, pnode + 1])
        edges.append([2 * n - 1, n - 1])
        edges.append([2 * n - 1, n])
        return Phylo(edge=edges, tip_label=labels, Nnode=n - 1)
    else:
        # Balanced tree
        edges = []
        next_internal = n + 1
        queue = [(1, n, next_internal)]
        next_internal += 1
        while queue:
            lo, hi, pnode = queue.pop(0)
            mid = (lo + hi) // 2
            if lo == mid:
                edges.append([pnode, lo])
            else:
                c1 = next_internal
                next_internal += 1
                edges.append([pnode, c1])
                queue.append((lo, mid, c1))
            if mid + 1 == hi:
                edges.append([pnode, hi])
            else:
                c2 = next_internal
                next_internal += 1
                edges.append([pnode, c2])
                queue.append((mid + 1, hi, c2))
        return Phylo(edge=edges, tip_label=labels, Nnode=next_internal - n - 1)


def estimate_dates(tree: Phylo, node_dates: Sequence[float], mu: float = 1.0, nsteps: int = 1) -> Dict[int, float]:
    """
    Estimate internal node dates from tip dates (compatible with ape::estimate.dates
    used in man/ggdensitree.Rd).
    """
    ntip = len(tree.tip_label)
    nnode_total = ntip + tree.Nnode
    dates: Dict[int, float] = {}
    for i, d in enumerate(node_dates, 1):
        dates[i] = float(d)
    post = reorder_phylo(tree, "postorder")
    children_of: Dict[int, List[int]] = {}
    for i in range(post.edge.shape[0]):
        p, ch = int(post.edge[i, 0]), int(post.edge[i, 1])
        children_of.setdefault(p, []).append(ch)
    for p, chs in children_of.items():
        ch_dates = [dates.get(c_node, 0.0) for c_node in chs]
        dates[p] = min(ch_dates) - 0.5 / max(mu, 1e-6)
    # Return 1-based indexable list/dict so dates[edge] works in Python/numpy
    arr = np.zeros(nnode_total + 1, dtype=float)
    for k, v in dates.items():
        if 0 <= k <= nnode_total:
            arr[k] = v
    return arr  # type: ignore[return-value]


# Polytomy & Clade extraction utilities
def di2multi(phy: Any, tol: float = 1e-8) -> Phylo:
    """
    Collapse branches shorter than tol into polytomies (multifurcations),
    compatible with ape::di2multi.
    """
    p = as_phylo(phy).copy()
    if p.edge_length is None:
        return p
    ntip = len(p.tip_label)
    edge = p.edge.copy()
    el = list(p.edge_length)

    # Collapse internal edges with length <= tol
    to_remove = []
    remap: Dict[int, int] = {}
    for i in range(edge.shape[0]):
        u, v = int(edge[i, 0]), int(edge[i, 1])
        if v > ntip and el[i] <= tol:
            # Resolve chained remappings
            target = remap.get(u, u)
            remap[v] = target
            to_remove.append(i)

    if not to_remove:
        return p

    def resolve(nd: int) -> int:
        while nd in remap:
            nd = remap[nd]
        return nd

    kept_edges: List[List[int]] = []
    kept_el: List[float] = []
    remove_set = set(to_remove)
    for i in range(edge.shape[0]):
        if i in remove_set:
            continue
        u = resolve(int(edge[i, 0]))
        v = int(edge[i, 1])
        kept_edges.append([u, v])
        kept_el.append(el[i])

    # Renumber internal nodes consecutively starting at ntip + 1
    internals = sorted({u for u, _ in kept_edges} | {v for _, v in kept_edges if v > ntip})
    id_map = {nd: nd for nd in range(1, ntip + 1)}
    for idx, old_id in enumerate(internals, ntip + 1):
        id_map[old_id] = idx

    final_edges = [[id_map[u], id_map[v]] for u, v in kept_edges]
    return Phylo(
        edge=final_edges,
        tip_label=list(p.tip_label),
        Nnode=len(internals),
        edge_length=kept_el,
        root_edge=p.root_edge,
        attrs=p.attrs,
    )


def as_polytomy(tree: Any, feature: str = "node.label", fun: Callable[[Any], bool] = lambda x: float(x) < 50) -> Phylo:
    """
    Collapse binary tree to polytomy by applying 'fun' to 'feature' (exact port of R/as.polytomy.R).
    """
    phy = as_phylo(tree).copy()
    ntip = len(phy.tip_label)
    if feature in ("branch.length", "edge.length"):
        if phy.edge_length is None:
            return phy
        to_collapse_edges = [i for i, val in enumerate(phy.edge_length) if not _is_scalar_na(val) and fun(val)]
    elif feature in ("node.label", "node_label"):
        if phy.node_label is None:
            return phy
        bad_nodes = set()
        for idx, val in enumerate(phy.node_label):
            if val is not None and not _is_scalar_na(val) and str(val) != "":
                try:
                    if fun(float(val)):
                        bad_nodes.add(ntip + 1 + idx)
                except Exception:
                    pass
        to_collapse_edges = [i for i in range(phy.edge.shape[0]) if int(phy.edge[i, 1]) in bad_nodes]
    else:
        return phy

    if not to_collapse_edges:
        return phy

    edge = phy.edge.copy()
    el = list(phy.edge_length) if phy.edge_length is not None else None
    collapse_set = set(to_collapse_edges)
    remap: Dict[int, int] = {}
    for idx in to_collapse_edges:
        p_node = int(edge[idx, 0])
        c_node = int(edge[idx, 1])
        if c_node <= ntip:
            continue
        remap[c_node] = p_node

    def _resolve_parent(nd: int) -> int:
        seen = set()
        cur = nd
        while cur in remap and cur not in seen:
            seen.add(cur)
            cur = remap[cur]
        return cur

    kept_edges: List[List[int]] = []
    kept_el_list: List[float] = []
    for i in range(edge.shape[0]):
        c_node = int(edge[i, 1])
        if i in collapse_set and c_node > ntip:
            continue
        p_node = _resolve_parent(int(edge[i, 0]))
        kept_edges.append([p_node, c_node])
        if el is not None:
            kept_el_list.append(el[i])

    kept_el = kept_el_list if el is not None else None
    internals = sorted({u for u, _ in kept_edges} | {v for _, v in kept_edges if v > ntip})
    id_map = {nd: nd for nd in range(1, ntip + 1)}
    for idx, old_id in enumerate(internals, ntip + 1):
        id_map[old_id] = idx
    final_edges = [[id_map[u], id_map[v]] for u, v in kept_edges]
    kept_nl = None
    if phy.node_label is not None:
        kept_nl = [
            phy.node_label[old_id - ntip - 1]
            if 0 <= (old_id - ntip - 1) < len(phy.node_label)
            else None
            for old_id in internals
        ]
    return Phylo(
        edge=final_edges,
        tip_label=list(phy.tip_label),
        Nnode=len(internals),
        edge_length=kept_el,
        node_label=kept_nl,
        root_edge=phy.root_edge,
        attrs=phy.attrs,
    )


def extract_clade(phy: Any, node: Union[int, str]) -> Phylo:
    """Extract subtree rooted at 'node' (compatible with ape::extract.clade)."""
    from .tidytree import nodeid, offspring
    p = as_phylo(phy)
    if isinstance(node, str):
        nid = nodeid(p, node)
    else:
        nid = int(node)
    ntip = len(p.tip_label)
    desc = offspring(p, nid, self_include=True)
    desc_set = set(int(x) for x in desc)
    tip_ids = sorted(x for x in desc_set if x <= ntip)
    int_ids = [nid] + sorted(x for x in desc_set if x > ntip and x != nid)

    id_map = {}
    for i, tid in enumerate(tip_ids, 1):
        id_map[tid] = i
    for i, iid in enumerate(int_ids, len(tip_ids) + 1):
        id_map[iid] = i

    new_edges = []
    new_el = []
    for i in range(p.edge.shape[0]):
        u, v = int(p.edge[i, 0]), int(p.edge[i, 1])
        if u in desc_set and v in desc_set:
            new_edges.append([id_map[u], id_map[v]])
            if p.edge_length is not None:
                new_el.append(float(p.edge_length[i]))

    new_nl = None
    if p.node_label is not None:
        new_nl = [
            p.node_label[iid - ntip - 1] if 0 <= (iid - ntip - 1) < len(p.node_label) else None
            for iid in int_ids
        ]

    return Phylo(
        edge=new_edges,
        tip_label=[p.tip_label[tid - 1] for tid in tip_ids],
        Nnode=len(int_ids),
        edge_length=new_el if p.edge_length is not None else None,
        node_label=new_nl,
    )


def system_file(*parts: str, package: Optional[str] = None) -> str:
    """
    Locate bundled example data files (compatible with R's system.file(..., package='treeio'/'ggtree')).
    """
    extdata_dir = os.path.join(os.path.dirname(__file__), "extdata")
    rel = os.path.join(*parts) if parts else ""
    basename = os.path.basename(rel)
    if basename == "sample.fasta":
        candidate_fas = os.path.join(extdata_dir, "pa.fas")
        if os.path.exists(candidate_fas):
            return candidate_fas
    candidate1 = os.path.join(extdata_dir, basename)
    if os.path.exists(candidate1):
        return candidate1
    candidate2 = os.path.join(extdata_dir, rel)
    if os.path.exists(candidate2):
        return candidate2
    # Check system R library if present
    if package:
        r_lib = f"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/{package}"
        candidate3 = os.path.join(r_lib, rel)
        if os.path.exists(candidate3):
            return candidate3
    return candidate1


def get_seqs(x: Any) -> Dict[str, str]:
    """Retrieve associated tip sequences from a TreeData object or FASTA path."""
    if isinstance(x, str):
        return read_fasta(x)
    if isinstance(x, TreeData):
        return dict(x.tip_seq or {})
    return {}


def subtrees(tree: Any) -> List[Phylo]:
    """Return list of all clades/subtrees of a phylogenetic tree (compatible with ape::subtrees)."""
    p = as_phylo(tree)
    ntip = len(p.tip_label)
    return [extract_clade(p, node_id) for node_id in range(ntip + 1, ntip + p.Nnode + 1)]

