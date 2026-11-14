"""
Port of tidytree and ggtree tree-traversal & grouping methods:
rootnode, isTip, nodeid, nodelab, child, parent, sibling, ancestor,
offspring, offspring_bulk, MRCA, groupOTU, groupClade, tree_subset,
td_filter, td_mutate, td_unnest.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Union
import numpy as np

from .dataframe import (
    ColumnSeries,
    DataFrame,
    Factor,
    TblTree,
    _is_scalar_na,
    c,
    factor,
    filter as df_filter,
    mutate as df_mutate,
    unnest as df_unnest,
)
from .tree import Phylo, TreeData, as_phylo, as_tibble, extract_clade, getNodeNum


def _is_ggtree_plot(x: Any) -> bool:
    return hasattr(x, "data") and isinstance(getattr(x, "data", None), DataFrame) and hasattr(x, "layers")


def taxa2node(data: DataFrame, taxa: Union[int, str]) -> int:
    """Resolve a taxon label or node ID to an integer node ID in a fortified DataFrame."""
    if isinstance(taxa, (int, np.integer)):
        return int(taxa)
    s = str(taxa)
    nodes = data["node"].tolist()
    labels = data["label"].tolist() if "label" in data else []
    for n, lbl in zip(nodes, labels):
        if lbl is not None and str(lbl) == s:
            return int(n)
    if s.isdigit():
        return int(s)
    raise ValueError(f"Taxon/node '{taxa}' not found in tree")


def nodeid(tree: Any, label: Union[str, Sequence[str]]) -> Union[int, List[Optional[int]]]:
    """
    Convert tip/internal node label(s) to 1-based integer node ID(s)
    (compatible with tidytree::nodeid and ggtree::nodeid.ggtree).
    """
    if _is_ggtree_plot(tree):
        return nodeid(tree.data, label)
    df = tree if isinstance(tree, DataFrame) else as_tibble(tree)
    labels_col = df["label"].tolist() if "label" in df else []
    nodes_col = [int(n) for n in df["node"].tolist()]
    lookup = {str(lbl): n for n, lbl in zip(nodes_col, labels_col) if not _is_scalar_na(lbl)}

    if isinstance(label, (str, int, np.integer)):
        if isinstance(label, (int, np.integer)):
            return int(label)
        return lookup.get(str(label), None)  # type: ignore[return-value]
    return [lookup.get(str(x), None) for x in label]


def nodelab(tree: Any, id: Union[int, Sequence[int]]) -> Union[Optional[str], List[Optional[str]]]:
    """
    Convert 1-based integer node ID(s) to node/tip label(s)
    (compatible with tidytree::nodelab and ggtree::nodelab.ggtree).
    """
    if _is_ggtree_plot(tree):
        return nodelab(tree.data, id)
    df = tree if isinstance(tree, DataFrame) else as_tibble(tree)
    labels_col = df["label"].tolist() if "label" in df else []
    nodes_col = [int(n) for n in df["node"].tolist()]
    lookup = {n: (str(lbl) if not _is_scalar_na(lbl) else None) for n, lbl in zip(nodes_col, labels_col)}

    if isinstance(id, (int, np.integer)):
        return lookup.get(int(id), None)
    return [lookup.get(int(x), None) for x in id]


def rootnode(x: Any) -> Any:
    """
    Find root node of a tree (returns 1-row TblTree for DataFrame, int node ID for Phylo/TreeData).
    """
    if _is_ggtree_plot(x):
        return rootnode(x.data)
    if isinstance(x, DataFrame):
        nodes = x["node"].tolist()
        parents = x["parent"].tolist()
        for i, (n, p) in enumerate(zip(nodes, parents)):
            if n == p:
                return x._slice_rows([i])
        # Fallback: parent not in node or first internal node
        node_set = set(nodes)
        for i, p in enumerate(parents):
            if p not in node_set:
                return x._slice_rows([i])
        return x._slice_rows([0])
    phy = as_phylo(x)
    edge = phy.edge
    parents = set(edge[:, 0].tolist())
    children = set(edge[:, 1].tolist())
    roots = list(parents - children)
    return int(roots[0]) if roots else len(phy.tip_label) + 1


def getRoot(x: Any) -> int:
    """Return integer root node ID."""
    r = rootnode(x)
    if isinstance(r, DataFrame):
        return int(r["node"][0])
    return int(r)


def isTip(x: Any, node: Optional[Union[int, Sequence[int]]] = None) -> Any:
    """Check whether node(s) are tips."""
    if _is_ggtree_plot(x):
        return isTip(x.data, node)
    df = x if isinstance(x, DataFrame) else as_tibble(x)
    parents = set(df["parent"].tolist())
    nodes = [int(n) for n in df["node"].tolist()]
    if "isTip" in df:
        tip_map = {int(n): bool(t) for n, t in zip(nodes, df["isTip"].tolist())}
    else:
        tip_map = {n: (n not in parents) for n in nodes}
    if node is None:
        return [tip_map[n] for n in nodes]
    if isinstance(node, (int, np.integer)):
        return tip_map.get(int(node), False)
    return [tip_map.get(int(n), False) for n in node]


def itself(data: DataFrame, node: Union[int, str]) -> DataFrame:
    """Return the 1-row DataFrame for 'node' (by integer ID or string label)."""
    nid = taxa2node(data, node)
    nodes = [int(n) for n in data["node"].tolist()]
    idxs = [i for i, n in enumerate(nodes) if n == nid]
    return data._slice_rows(idxs)


def child(data: Any, node: Union[int, str]) -> Any:
    """
    Get direct children of 'node'.
    Returns DataFrame if data is DataFrame/GGTree, or list of int node IDs if Phylo/TreeData.
    """
    if _is_ggtree_plot(data):
        return child(data.data, node)
    if isinstance(data, DataFrame):
        nid = taxa2node(data, node)
        nodes = [int(n) for n in data["node"].tolist()]
        parents = [int(p) for p in data["parent"].tolist()]
        idxs = [i for i, (p, n) in enumerate(zip(parents, nodes)) if p == nid and p != n]
        return data._slice_rows(idxs)
    phy = as_phylo(data)
    nid = int(nodeid(phy, node)) if isinstance(node, str) else int(node)
    edge = phy.edge
    return [int(edge[i, 1]) for i in range(edge.shape[0]) if int(edge[i, 0]) == nid and int(edge[i, 1]) != nid]


def parent(data: Any, node: Union[int, str]) -> Any:
    """
    Get parent of 'node'.
    Returns DataFrame if data is DataFrame/GGTree, or int node ID (0 if root) if Phylo/TreeData.
    """
    if _is_ggtree_plot(data):
        return parent(data.data, node)
    if isinstance(data, DataFrame):
        nd_df = itself(data, node)
        if len(nd_df) == 0:
            return data._slice_rows([])
        nid = int(nd_df["node"][0])
        pid = int(nd_df["parent"][0])
        if pid == nid:
            return data._slice_rows([])
        nodes = [int(n) for n in data["node"].tolist()]
        idxs = [i for i, n in enumerate(nodes) if n == pid]
        return data._slice_rows(idxs)
    phy = as_phylo(data)
    nid = int(nodeid(phy, node)) if isinstance(node, str) else int(node)
    edge = phy.edge
    for i in range(edge.shape[0]):
        if int(edge[i, 1]) == nid:
            return int(edge[i, 0])
    return 0


def sibling(data: Any, node: Union[int, str]) -> Any:
    """Get sibling(s) of 'node'."""
    if _is_ggtree_plot(data):
        return sibling(data.data, node)
    if isinstance(data, DataFrame):
        nid = taxa2node(data, node)
        p_df = parent(data, nid)
        if len(p_df) == 0:
            return p_df
        pid = int(p_df["node"][0])
        ch_df = child(data, pid)
        nodes = [int(n) for n in ch_df["node"].tolist()]
        idxs = [i for i, n in enumerate(nodes) if n != nid]
        return ch_df._slice_rows(idxs)
    phy = as_phylo(data)
    nid = int(nodeid(phy, node)) if isinstance(node, str) else int(node)
    pid = parent(phy, nid)
    if pid == 0:
        return []
    return [c_id for c_id in child(phy, pid) if c_id != nid]


def ancestor(data: Any, node: Union[int, str]) -> Any:
    """
    Get all ancestors of 'node' up to the root.
    Returns DataFrame if data is DataFrame, or list of int node IDs if Phylo/TreeData.
    """
    if _is_ggtree_plot(data):
        return ancestor(data.data, node)
    if isinstance(data, DataFrame):
        nd_df = itself(data, node)
        if len(nd_df) == 0:
            return data._slice_rows([])
        nid = int(nd_df["node"][0])
        pid = int(nd_df["parent"][0])
        if nid == pid:
            return data._slice_rows([])
        nodes = [int(n) for n in data["node"].tolist()]
        parents = [int(p) for p in data["parent"].tolist()]
        p_map = {n: p for n, p in zip(nodes, parents)}
        anc_ids = []
        cur = pid
        visited = set()
        while cur not in visited:
            visited.add(cur)
            anc_ids.append(cur)
            nxt = p_map.get(cur, cur)
            if nxt == cur:
                break
            cur = nxt
        anc_set = set(anc_ids)
        idxs = [i for i, n in enumerate(nodes) if n in anc_set]
        return data._slice_rows(idxs)

    df = as_tibble(data)
    anc_df = ancestor(df, node)
    return [int(n) for n in anc_df["node"].tolist()]


def offspring(
    data: Any,
    node: Union[int, str, Sequence[Union[int, str]]],
    tiponly: bool = False,
    self_include: bool = False,
    type: str = "all",
    **kwargs: Any,
) -> Any:
    """
    Get all descendant (offspring) nodes of 'node'
    (exact port of tidytree::offspring.tbl_tree, offspring.phylo, and ggtree::offspring.ggtree).
    """
    if _is_ggtree_plot(data):
        return offspring(data.data, node, tiponly=tiponly, self_include=self_include, type=type, **kwargs)

    if isinstance(node, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(node, (str, bytes)):
        node_list = list(node)
        if len(node_list) == 1:
            return offspring(data, node_list[0], tiponly=tiponly, self_include=self_include, type=type, **kwargs)
        return {
            n: offspring(data, n, tiponly=tiponly, self_include=self_include, type=type, **kwargs)
            for n in node_list
        }

    if isinstance(data, DataFrame):
        nid = taxa2node(data, node)  # type: ignore[arg-type]
        nodes = [int(n) for n in data["node"].tolist()]
        parents = [int(p) for p in data["parent"].tolist()]
        rn = getRoot(data)

        kids: Dict[int, List[int]] = {}
        for p, ch in zip(parents, nodes):
            if p != ch and ch != rn:
                kids.setdefault(p, []).append(ch)

        direct = kids.get(nid, [])
        if not direct:
            if self_include:
                idxs = [i for i, n in enumerate(nodes) if n == nid]
                return data._slice_rows(idxs)
            return data._slice_rows([])

        desc_ids = list(direct)
        i = 0
        while i < len(desc_ids):
            cur = desc_ids[i]
            if cur in kids:
                desc_ids.extend(kids[cur])
            i += 1

        if self_include:
            desc_ids = [nid] + desc_ids

        desc_set = set(desc_ids)
        if tiponly or type in ("tips", "external"):
            if "isTip" in data:
                is_tip_col = [bool(t) for t in data["isTip"].tolist()]
                idxs = [i for i, (n, t) in enumerate(zip(nodes, is_tip_col)) if n in desc_set and t]
            else:
                idxs = [i for i, n in enumerate(nodes) if n in desc_set and n < rn]
        elif type == "internal":
            idxs = [i for i, n in enumerate(nodes) if n in desc_set and n >= rn]
        elif type == "children":
            direct_set = set(direct) | ({nid} if self_include else set())
            idxs = [i for i, n in enumerate(nodes) if n in direct_set]
        else:
            idxs = [i for i, n in enumerate(nodes) if n in desc_set]
        return data._slice_rows(idxs)

    # Phylo / TreeData returns list of int node IDs
    df = as_tibble(data)
    sub = offspring(df, node, tiponly=tiponly, self_include=self_include, type=type, **kwargs)
    return [int(n) for n in sub["node"].tolist()]


def offspring_bulk(data: DataFrame, nodes: Sequence[Union[int, str]], self_include: bool = False) -> List[List[int]]:
    """
    Bulk offspring traversal returning a list of integer node ID vectors for each node in 'nodes'
    (exact port of ggtree's internal offspring_bulk in R/utilities.R).
    """
    node_col = [int(n) for n in data["node"].tolist()]
    parent_col = [int(p) for p in data["parent"].tolist()]
    rn = getRoot(data)
    kids: Dict[int, List[int]] = {}
    for p, ch in zip(parent_col, node_col):
        if p != ch and ch != rn:
            kids.setdefault(p, []).append(ch)

    result: List[List[int]] = []
    for raw_nd in nodes:
        nid = taxa2node(data, raw_nd)
        stack = list(kids.get(nid, []))
        desc = [nid] if self_include else []
        while stack:
            cur = stack.pop(0)
            desc.append(cur)
            if cur in kids:
                stack.extend(kids[cur])
        result.append(desc)
    return result


def _mrca_pair_df(data: DataFrame, node1: Union[int, str], node2: Union[int, str]) -> DataFrame:
    nid1 = taxa2node(data, node1)
    nid2 = taxa2node(data, node2)
    if nid1 == nid2:
        return itself(data, nid1)
    anc1 = ancestor(data, nid1)
    if len(anc1) == 0:
        return itself(data, nid1)
    anc1_ids = set(int(n) for n in anc1["node"].tolist()) | {nid1}
    if nid2 in anc1_ids:
        return itself(data, nid2)
    p = parent(data, nid2)
    while len(p) > 0 and int(p["node"][0]) not in anc1_ids:
        p = parent(data, int(p["node"][0]))
    return p


class MRCARow(TblTree):
    """
    1-row TblTree returned by MRCA(tbl_tree, ...) that also behaves like an int
    when compared or cast to int (so both `MRCA(df, a, b)$node` and `int(MRCA(df, a, b))` work!).
    """

    def __int__(self) -> int:
        return int(self["node"][0]) if len(self) > 0 else 0

    def __index__(self) -> int:
        return self.__int__()

    @property
    def node(self) -> Any:
        # In R: MRCA(df, a, b)$node returns scalar int when 1 row!
        if len(self) == 1:
            return int(self._cols["node"][0])
        return ColumnSeries(self._cols.get("node", []), name="node")

    def __eq__(self, other: Any) -> Any:  # type: ignore[override]
        if isinstance(other, (int, np.integer)):
            return self.__int__() == int(other)
        return super().__eq__(other)


def MRCA(data: Any, node1: Any = None, node2: Any = None, *more_nodes: Any, **kwargs: Any) -> Any:
    """
    Find Most Recent Common Ancestor of 2 or more nodes/taxa.
    - For Phylo / TreeData / GGTree plot: returns integer node ID (matching tidytree::MRCA.phylo and ggtree::MRCA.ggtree).
    - For DataFrame / TblTree: returns MRCARow (1-row TblTree where .node and int() both return the node ID!).
    """
    if node1 is None and "tip" in kwargs:
        node1 = kwargs.pop("tip")
    if node1 is None and ".node1" in kwargs:
        node1 = kwargs.pop(".node1")
    if node2 is None and ".node2" in kwargs:
        node2 = kwargs.pop(".node2")
    if _is_ggtree_plot(data):
        res_df = MRCA(data.data, node1, node2, *more_nodes)
        return int(res_df["node"][0]) if len(res_df) > 0 else None

    if not isinstance(data, DataFrame):
        df = as_tibble(data)
        res_df = MRCA(df, node1, node2, *more_nodes)
        return int(res_df["node"][0]) if len(res_df) > 0 else None

    # Gather all target nodes
    targets: List[Union[int, str]] = []
    for arg in (node1, node2, *more_nodes):
        if arg is None:
            continue
        if isinstance(arg, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(arg, (str, bytes)):
            targets.extend(list(arg))
        else:
            targets.append(arg)

    if len(targets) == 0:
        raise ValueError("MRCA requires at least one node")
    if len(targets) == 1:
        row = itself(data, targets[0])
        return MRCARow(row, attrs=data.attrs)

    cur_df = _mrca_pair_df(data, targets[0], targets[1])
    for nxt in targets[2:]:
        if len(cur_df) == 0:
            break
        cur_df = _mrca_pair_df(data, int(cur_df["node"][0]), nxt)
    return MRCARow(cur_df, attrs=data.attrs)


# groupOTU and groupClade (exact ports of tidytree::groupOTU and tidytree::groupClade + ggtree methods)
def _groupOTU_tbl_tree_item(
    data: DataFrame,
    focus_nodes: Sequence[Union[int, str]],
    focus_label: Optional[Any] = None,
    group_name: str = "group",
    overlap: str = "overwrite",
    connect: bool = False,
) -> DataFrame:
    out = data.copy()
    nodes_col = [int(n) for n in out["node"].tolist()]
    labels_col = out["label"].tolist() if "label" in out else [None] * len(nodes_col)

    focus_ids: List[int] = []
    for item in focus_nodes:
        if isinstance(item, (int, np.integer)):
            focus_ids.append(int(item))
        else:
            s = str(item)
            for n, lbl in zip(nodes_col, labels_col):
                if lbl is not None and str(lbl) == s:
                    focus_ids.append(n)

    n = len(out)
    if group_name not in out:
        foc: List[Any] = [0] * n
    else:
        foc = list(out[group_name].tolist())

    # Determine next numeric group index g
    numeric_vals = [0]
    for v in foc:
        if not _is_scalar_na(v):
            try:
                numeric_vals.append(int(float(v)))
            except Exception:
                pass
    g = max(numeric_vals) + 1
    if focus_label is None:
        focus_label = g

    if not focus_ids:
        out[group_name] = foc
        return out

    if len(focus_ids) == 1:
        hit = list(focus_ids)
    else:
        anc_lists = [sorted(int(x) for x in ancestor(out, fid)["node"].tolist()) for fid in focus_ids]
        min_len = min((len(a) for a in anc_lists), default=0)
        i = 1  # 0-based index corresponding to R's i <- 2L
        while i < min_len:
            vals_at_i = {a[i] for a in anc_lists}
            if len(vals_at_i) != 1:
                break
            i += 1
        # Drop 0..(i-1) except keep anc_lists[0][i-1] (the MRCA of focus)
        common_mrca = [anc_lists[0][i - 1]] if i - 1 >= 0 and len(anc_lists[0]) >= i else []
        remaining_anc = [x for a in anc_lists for x in a[i:]]
        hit = list(dict.fromkeys(common_mrca + remaining_anc + focus_ids))

    # Map node ID -> row index in out
    id_to_row = {nd: idx for idx, nd in enumerate(nodes_col)}

    if overlap == "origin":
        sn = [nid for nid in hit if nid in id_to_row and (_is_scalar_na(foc[id_to_row[nid]]) or foc[id_to_row[nid]] == 0 or foc[id_to_row[nid]] == "0")]
    elif overlap == "abandon":
        sn = []
        for nid in hit:
            if nid not in id_to_row:
                continue
            ridx = id_to_row[nid]
            if not _is_scalar_na(foc[ridx]) and foc[ridx] != 0 and foc[ridx] != "0":
                foc[ridx] = None
            else:
                sn.append(nid)
    else:
        sn = [nid for nid in hit if nid in id_to_row]

    if sn and connect:
        parents_col = [int(p) for p in out["parent"].tolist()]
        deg_counts: Dict[int, int] = {}
        for nid in sn:
            p = parents_col[id_to_row[nid]]
            deg_counts[p] = deg_counts.get(p, 0) + 1
        if sum(1 for cnt in deg_counts.values() if cnt > 1) == 1:
            sn = [nid for nid in focus_ids if nid in id_to_row]

    for nid in sn:
        foc[id_to_row[nid]] = focus_label

    out[group_name] = foc
    return out


def groupOTU(
    data: Any,
    _node: Any,
    group_name: str = "group",
    overlap: str = "overwrite",
    connect: bool = False,
    **kwargs: Any,
) -> Any:
    """
    Group OTUs by tracing back to their most recent common ancestor
    (compatible with tidytree::groupOTU and ggtree::groupOTU.ggtree).
    Supports Phylo, TreeData, DataFrame/TblTree, and GGTree plot objects.
    """
    if _is_ggtree_plot(data):
        out_plot = data.copy()
        out_plot.data = groupOTU(out_plot.data, _node, group_name=group_name, overlap=overlap, connect=connect, **kwargs)
        if group_name not in out_plot.mapping:
            from .aes import aes
            out_plot = out_plot + aes(color=group_name)
        return out_plot

    if isinstance(data, TreeData):
        out_td = data.copy()
        out_td.phylo = groupOTU(out_td.phylo, _node, group_name=group_name, overlap=overlap, connect=connect, **kwargs)
        return out_td

    if isinstance(data, Phylo):
        phy = data.copy()
        df = as_tibble(phy)
        df_grouped = groupOTU(df, _node, group_name=group_name, overlap=overlap, connect=connect, **kwargs)
        df_sorted = df_grouped.sort_values("node")
        vals = [0 if _is_scalar_na(v) else v for v in df_sorted[group_name].tolist()]
        phy.attrs[group_name] = factor(vals)
        return phy

    if isinstance(data, DataFrame):
        out = data.copy()
        out[group_name] = None  # Delete existing column if present
        if isinstance(_node, dict):
            for k, v in _node.items():
                v_list = list(v) if isinstance(v, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(v, (str, bytes)) else [v]
                out = _groupOTU_tbl_tree_item(out, v_list, focus_label=k, group_name=group_name, overlap=overlap, connect=connect)
        elif isinstance(_node, (list, tuple)) and len(_node) > 0 and isinstance(_node[0], (list, tuple, np.ndarray)):
            for idx, v in enumerate(_node, 1):
                out = _groupOTU_tbl_tree_item(out, list(v), focus_label=idx, group_name=group_name, overlap=overlap, connect=connect)
        else:
            v_list = list(_node) if isinstance(_node, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(_node, (str, bytes)) else [_node]
            out = _groupOTU_tbl_tree_item(out, v_list, focus_label=None, group_name=group_name, overlap=overlap, connect=connect)

        rn = getRoot(out)
        nodes_col = [int(n) for n in out["node"].tolist()]
        if rn in nodes_col:
            rn_idx = nodes_col.index(rn)
            col_vals = out[group_name].tolist()
            rn_val = col_vals[rn_idx]
            if sum(1 for v in col_vals if v == rn_val) == 1:
                col_vals[rn_idx] = 0
                out[group_name] = col_vals
        out[group_name] = factor(out[group_name].tolist())
        return out

    raise TypeError(f"Unsupported type for groupOTU: {type(data)}")


def groupClade(
    data: Any,
    _node: Any,
    group_name: str = "group",
    overlap: str = "overwrite",
    **kwargs: Any,
) -> Any:
    """
    Group clade(s) rooted at '_node'
    (compatible with tidytree::groupClade and ggtree::groupClade.ggtree).
    Supports Phylo, TreeData, DataFrame/TblTree, and GGTree plot objects.
    """
    if _is_ggtree_plot(data):
        out_plot = data.copy()
        out_plot.data = groupClade(out_plot.data, _node, group_name=group_name, overlap=overlap, **kwargs)
        if group_name not in out_plot.mapping:
            from .aes import aes
            out_plot = out_plot + aes(color=group_name)
        return out_plot

    if isinstance(data, TreeData):
        out_td = data.copy()
        out_td.phylo = groupClade(out_td.phylo, _node, group_name=group_name, overlap=overlap, **kwargs)
        return out_td

    if isinstance(data, Phylo):
        phy = data.copy()
        df = as_tibble(phy)
        df_grouped = groupClade(df, _node, group_name=group_name, overlap=overlap, **kwargs)
        df_sorted = df_grouped.sort_values("node")
        phy.attrs[group_name] = factor(df_sorted[group_name].tolist())
        return phy

    if isinstance(data, DataFrame):
        out = data.copy()
        n = len(out)
        nodes_col = [int(nd) for nd in out["node"].tolist()]
        id_to_row = {nd: idx for idx, nd in enumerate(nodes_col)}
        foc: List[Any] = [0] * n

        if isinstance(_node, dict):
            node_items = list(_node.items())
        elif isinstance(_node, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(_node, (str, bytes)):
            node_items = [(None, x) for x in _node]
        else:
            node_items = [(None, _node)]

        if len(node_items) == 1 and node_items[0][0] is None:
            nid = taxa2node(out, node_items[0][1])
            hit = [nid] + [int(x) for x in offspring(out, nid)["node"].tolist()]
            for h in hit:
                if h in id_to_row:
                    foc[id_to_row[h]] = 1
            out[group_name] = factor(foc)
            return out

        for name_k, raw_nd in node_items:
            nid = taxa2node(out, raw_nd)
            hit = [nid] + [int(x) for x in offspring(out, nid)["node"].tolist()]
            if overlap == "origin":
                sn = [h for h in hit if h in id_to_row and (_is_scalar_na(foc[id_to_row[h]]) or foc[id_to_row[h]] == 0)]
            elif overlap == "abandon":
                sn = []
                for h in hit:
                    if h not in id_to_row:
                        continue
                    ridx = id_to_row[h]
                    if not _is_scalar_na(foc[ridx]) and foc[ridx] != 0:
                        foc[ridx] = None
                    else:
                        sn.append(h)
            else:
                sn = [h for h in hit if h in id_to_row]

            if sn:
                if name_k is None:
                    num_vals = [0]
                    for v in foc:
                        if not _is_scalar_na(v):
                            try:
                                num_vals.append(int(float(v)))
                            except Exception:
                                pass
                    lbl_val: Any = max(num_vals) + 1
                else:
                    lbl_val = name_k
                for h in sn:
                    foc[id_to_row[h]] = lbl_val

        out[group_name] = factor(foc)
        return out

    raise TypeError(f"Unsupported type for groupClade: {type(data)}")


def tree_subset(tree: Any, node: Union[int, str], levels_back: int = 5, group_node: bool = True) -> TreeData:
    """Extract a local neighborhood subtree around 'node' going back 'levels_back' ancestors."""
    df = as_tibble(tree)
    nid = taxa2node(df, node)
    cur = nid
    for _ in range(int(levels_back)):
        p_df = parent(df, cur)
        if len(p_df) == 0:
            break
        cur = int(p_df["node"][0])
    sub_phy = extract_clade(tree, cur)
    td = TreeData(phylo=sub_phy)
    if group_node:
        focus_lbl = nodelab(df, nid)
        if focus_lbl:
            td = groupOTU(td, [focus_lbl])
    return td


# ggfun data selector helpers: td_filter, td_mutate, td_unnest
def td_filter(*conditions: Any, **kwargs: Any) -> Callable[[DataFrame], DataFrame]:
    """Return a function that filters a tree DataFrame (compatible with ggfun::td_filter)."""
    from .aes import _capture_caller_env
    captured_env = _capture_caller_env()

    def _fn(d: DataFrame) -> DataFrame:
        return df_filter(d, *conditions, _extra_env=captured_env, **kwargs)
    return _fn


def td_mutate(**kwargs: Any) -> Callable[[DataFrame], DataFrame]:
    """Return a function that mutates a tree DataFrame (compatible with ggfun::td_mutate)."""
    from .aes import _capture_caller_env
    captured_env = _capture_caller_env()

    def _fn(d: DataFrame) -> DataFrame:
        return df_mutate(d, _extra_env=captured_env, **kwargs)
    return _fn


def td_unnest(cols: Union[str, Sequence[str]], **kwargs: Any) -> Callable[[DataFrame], DataFrame]:
    """Return a function that unnests list-columns in a tree DataFrame (compatible with ggfun::td_unnest)."""
    def _fn(d: DataFrame) -> DataFrame:
        return df_unnest(d, cols)
    return _fn
