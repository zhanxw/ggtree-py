"""
Exact port of ggtree / treelayout tree coordinate & geometry algorithms:
fortify, tree_layout, all 15+ layouts (rectangular, slanted, circular, fan,
inward_circular, radial, unrooted, equal_angle, daylight, ape, tree_and_leaf,
dendrogram, ellipse, roundrect, tidy, external callables), clade_geometry,
calculate_angle, calculate_branch_mid, Date2decimal, decimal2Date,
get_layout, tip_y_alignment, etc.
"""

from __future__ import annotations

import datetime
import math
import warnings
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union
import numpy as np

from .dataframe import (
    ColumnSeries,
    DataFrame,
    TblTree,
    _is_scalar_na,
    bind_rows,
    ensure_dataframe,
    full_join,
    getOption,
    is_dataframe_like,
)
from .tidytree import (
    ancestor,
    child,
    getRoot,
    offspring,
    offspring_bulk,
    parent,
    rootnode,
    sibling,
    taxa2node,
)
from .tree import (
    Hclust,
    MultiPhylo,
    Nnode,
    Nnode2,
    Ntip,
    Phylo,
    TreeData,
    TreeDataList,
    as_phylo,
    as_tibble,
    getNodeNum,
    ladderize as ladderize_phylo,
    reorder_phylo,
)


# Date utilities
def _parse_date(d: Any) -> datetime.date:
    if isinstance(d, datetime.datetime):
        return d.date()
    if isinstance(d, datetime.date):
        return d
    s = str(d).strip()
    parts = s.split("-")
    return datetime.date(int(parts[0]), int(parts[1]), int(parts[2]))


class RDate(datetime.date):
    """Date object supporting string ISO comparison (R Date interoperability)."""

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.isoformat() == other
        return super().__eq__(other)

    def __ne__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.isoformat() != other
        return super().__ne__(other)

    __hash__ = datetime.date.__hash__


def _is_leap_year(year: int) -> bool:
    return (year % 4 == 0) and (year % 100 != 0 or year % 400 == 0)


def Date2decimal(x: Any) -> Any:
    """Convert Date or 'YYYY-MM-DD' string(s) to decimal year (exact port of ggtree::Date2decimal)."""
    if isinstance(x, (int, float, np.number)):
        return float(x)
    if isinstance(x, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(x, (str, bytes)):
        return [Date2decimal(item) for item in x]
    dt = _parse_date(x)
    year = dt.year
    jan1 = datetime.date(year, 1, 1)
    days = (dt - jan1).days
    ndate = 366.0 if _is_leap_year(year) else 365.0
    return float(year) + float(days) / ndate


def decimal2Date(x: Any) -> Any:
    """Convert decimal year(s) back to RDate / datetime.date (exact port of ggtree::decimal2Date)."""
    if isinstance(x, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(x, (str, bytes)):
        return [decimal2Date(item) for item in x]
    val = float(x)
    yr = int(math.floor(val))
    frac = val - yr
    jan1 = datetime.date(yr, 1, 1)
    ndate = 366.0 if _is_leap_year(yr) else 365.0
    days = int(math.floor(frac * ndate + 1e-9))
    res = jan1 + datetime.timedelta(days=days)
    return RDate(res.year, res.month, res.day)


# Rectangular / Slanted / Cladogram X and Y coordinate calculators
def getXcoord2(
    x: Sequence[Optional[float]],
    root: int,
    parent_arr: Sequence[int],
    child_arr: Sequence[int],
    len_arr: Sequence[float],
    start: float = 0.0,
    rev: bool = False,
) -> List[float]:
    """Exact port of ggtree:::getXcoord2 (1-based node indexing)."""
    n = len(x)
    out: List[Optional[float]] = [None] * n
    out[root - 1] = float(start)
    direction = -1.0 if rev else 1.0
    ignore_neg = bool(getOption("ignore.negative.edge", False))

    if any(l < 0 for l in len_arr if not _is_scalar_na(l)) and not ignore_neg:
        warnings.warn("The tree contained negative edge length(s).")

    parent_to_edges: Dict[int, List[Tuple[int, float]]] = {}
    for p, ch, l in zip(parent_arr, child_arr, len_arr):
        if p != ch:
            parent_to_edges.setdefault(int(p), []).append((int(ch), float(l) if not _is_scalar_na(l) else 0.0))

    current_nodes = [int(root)]
    visited = {int(root)}
    while current_nodes:
        new_nodes: List[int] = []
        for pnode in current_nodes:
            px = out[pnode - 1] or 0.0
            for ch, edge_l in parent_to_edges.get(pnode, []):
                if ch in visited:
                    continue
                visited.add(ch)
                if ignore_neg:
                    step = abs(edge_l) * direction
                else:
                    step = edge_l * direction
                out[ch - 1] = px + step
                new_nodes.append(ch)
        current_nodes = new_nodes

    return [float(v) if v is not None else 0.0 for v in out]


def getXcoord(tr: Phylo) -> List[float]:
    """Exact port of ggtree:::getXcoord."""
    edge = tr.edge
    parent_arr = edge[:, 0].tolist()
    child_arr = edge[:, 1].tolist()
    root = getRoot(tr)
    len_arr = tr.edge_length if tr.edge_length is not None else [1.0] * edge.shape[0]
    N = getNodeNum(tr)
    x = [0.0] * N
    return getXcoord2(x, root, parent_arr, child_arr, len_arr)


def getXcoord_no_length(tr: Phylo) -> List[float]:
    """Exact port of ggtree:::getXcoord_no_length (cladogram x coordinates)."""
    edge = tr.edge
    parent_arr = [int(v) for v in edge[:, 0].tolist()]
    child_arr = [int(v) for v in edge[:, 1].tolist()]
    N = getNodeNum(tr)
    ntip = Ntip(tr)

    x: List[Optional[float]] = [0.0 if i < ntip else None for i in range(N)]
    child_list: Dict[int, List[int]] = {}
    child_to_parent: Dict[int, int] = {}
    for p, ch in zip(parent_arr, child_arr):
        if p != ch:
            child_list.setdefault(p, []).append(ch)
            child_to_parent[ch] = p

    current_nodes = set(range(1, ntip + 1))
    while any(v is None for v in x):
        p_candidates = sorted({child_to_parent[c] for c in current_nodes if c in child_to_parent})
        new_nodes = [p for p in p_candidates if all(ch in current_nodes for ch in child_list.get(p, []))]
        if not new_nodes:
            break
        exclude = set()
        for j in new_nodes:
            chs = child_list[j]
            x[j - 1] = min(float(x[ch - 1]) for ch in chs if x[ch - 1] is not None) - 1.0
            exclude.update(chs)
        current_nodes = (current_nodes - exclude) | set(new_nodes)

    valid = [float(v) for v in x if v is not None]
    min_x = min(valid) if valid else 0.0
    return [(float(v) if v is not None else min_x) - min_x for v in x]


def _convert_tips2ancestors_sbp(tree: Phylo, include_root: bool = True) -> Tuple[np.ndarray, List[int]]:
    """Exact port of ggtree:::.convert_tips2ancestors_sbp for slanted cladograms."""
    ntip = Ntip(tree)
    N = getNodeNum(tree)
    root = getRoot(tree)
    all_nodes = list(range(1, N + 1))
    if not include_root and root in all_nodes:
        all_nodes.remove(root)

    df = as_tibble(tree)
    sbp = np.zeros((ntip, len(all_nodes)), dtype=float)
    col_index = {nid: idx for idx, nid in enumerate(all_nodes)}
    for tip_id in range(1, ntip + 1):
        anc = [int(x) for x in ancestor(df, tip_id)["node"].tolist()] + [tip_id]
        for a in anc:
            if a in col_index:
                sbp[tip_id - 1, col_index[a]] = 1.0
    return sbp, all_nodes


def getXcoord_no_length_slanted(sbp: np.ndarray, all_nodes: List[int]) -> List[float]:
    """Exact port of ggtree:::getXcoord_no_length_slanted."""
    col_sums = -np.sum(sbp, axis=0)
    max_abs = float(np.max(np.abs(col_sums))) if col_sums.size > 0 else 0.0
    x_by_node = {nid: float(col_sums[idx] + max_abs) for idx, nid in enumerate(all_nodes)}
    return [x_by_node[nid] for nid in sorted(all_nodes)]


def getYcoord_no_length_slanted(sbp: np.ndarray, all_nodes: List[int]) -> List[float]:
    """Exact port of ggtree:::getYcoord_no_length_slanted."""
    nrow = sbp.shape[0]
    row_weights = np.arange(1, nrow + 1, dtype=float).reshape(-1, 1)
    weighted = row_weights * sbp
    y_by_node: Dict[int, float] = {}
    for idx, nid in enumerate(all_nodes):
        col = weighted[:, idx]
        nz = col[col != 0]
        y_by_node[nid] = float(np.mean(nz)) if nz.size > 0 else 0.0
    return [y_by_node[nid] for nid in sorted(all_nodes)]


def getYcoord(tr: Phylo, step: float = 1.0, tip_order: Optional[Sequence[str]] = None) -> List[float]:
    """Exact port of ggtree:::getYcoord."""
    ntip = len(tr.tip_label)
    N = getNodeNum(tr)
    edge = tr.edge
    parent_arr = [int(v) for v in edge[:, 0].tolist()]
    child_arr = [int(v) for v in edge[:, 1].tolist()]

    child_list: Dict[int, List[int]] = {}
    pvec: Dict[int, int] = {}
    for p, ch in zip(parent_arr, child_arr):
        if p != ch:
            child_list.setdefault(p, []).append(ch)
            pvec[ch] = p

    y: List[Optional[float]] = [None] * N
    if tip_order is None:
        tip_in_edge_order = [ch for ch in child_arr if ch <= ntip]
        for idx, tid in enumerate(tip_in_edge_order, 1):
            y[tid - 1] = float(idx) * step
    else:
        order_map = {str(lbl): idx for idx, lbl in enumerate(tip_order, 1)}
        for tid in range(1, ntip + 1):
            lbl = tr.tip_label[tid - 1]
            y[tid - 1] = float(order_map.get(lbl, tid)) * step

    current_nodes = set(range(1, ntip + 1))
    while any(v is None for v in y):
        p_nodes = []
        for c_id in current_nodes:
            if c_id in pvec and pvec[c_id] not in p_nodes:
                p_nodes.append(pvec[c_id])
        new_nodes = [p for p in p_nodes if all(ch in current_nodes for ch in child_list.get(p, []))]
        if not new_nodes:
            break
        consumed = set()
        for p in new_nodes:
            chs = child_list[p]
            vals = [float(y[ch - 1]) for ch in chs if y[ch - 1] is not None]
            y[p - 1] = sum(vals) / len(vals) if vals else 0.0
            consumed.update(chs)
        current_nodes = (current_nodes - consumed) | set(new_nodes)

    return [float(v) if v is not None else 0.0 for v in y]


def re_assign_ycoord_df(df: DataFrame, currentNode: Sequence[int]) -> DataFrame:
    """Exact port of ggtree:::re_assign_ycoord_df."""
    out = df.copy()
    nodes = [int(n) for n in out["node"].tolist()]
    parents = [int(p) for p in out["parent"].tolist()]
    y_vals = list(out["y"].tolist())
    id_to_idx = {n: i for i, n in enumerate(nodes)}

    child_map: Dict[int, List[int]] = {}
    for p, n in zip(parents, nodes):
        if p != n:
            child_map.setdefault(p, []).append(n)

    cur_set = [int(c) for c in currentNode]
    while any(_is_scalar_na(v) for v in y_vals):
        p_nodes: List[int] = []
        for c_id in cur_set:
            if c_id in id_to_idx:
                p = parents[id_to_idx[c_id]]
                if p not in p_nodes:
                    p_nodes.append(p)
        new_nodes = [p for p in p_nodes if all(ch in cur_set for ch in child_map.get(p, []))]
        if not new_nodes:
            break
        traced = set()
        for p in new_nodes:
            chs = child_map.get(p, [])
            ch_ys = [float(y_vals[id_to_idx[ch]]) for ch in chs if not _is_scalar_na(y_vals[id_to_idx[ch]])]
            if p in id_to_idx and ch_ys:
                y_vals[id_to_idx[p]] = sum(ch_ys) / len(ch_ys)
            traced.update(chs)
        cur_set = [c for c in cur_set if c not in traced] + new_nodes

    out["y"] = y_vals
    return out


def reassign_y_from_node_to_root(df: DataFrame, node: int) -> DataFrame:
    """Exact port of ggtree:::reassign_y_from_node_to_root."""
    out = df.copy()
    nodes = [int(n) for n in out["node"].tolist()]
    parents = [int(p) for p in out["parent"].tolist()]
    y_vals: List[Optional[float]] = [None if _is_scalar_na(v) else float(v) for v in out["y"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes)}

    root_id = getRoot(out)
    cur = int(node)
    if cur not in id_to_idx:
        return out
    pp = parents[id_to_idx[cur]]
    visited = set()
    while pp != root_id and pp not in visited:
        visited.add(pp)
        ch_df = child(out, pp)
        ch_nodes = [int(n) for n in ch_df["node"].tolist()]
        ch_ys = [y_vals[id_to_idx[c]] for c in ch_nodes if c in id_to_idx and y_vals[id_to_idx[c]] is not None]
        if ch_ys and pp in id_to_idx:
            y_vals[id_to_idx[pp]] = sum(ch_ys) / len(ch_ys)  # type: ignore[arg-type]
        pp = parents[id_to_idx[pp]]

    root_chs = [int(n) for n in child(out, root_id)["node"].tolist() if int(n) != root_id]
    root_ys = [y_vals[id_to_idx[c]] for c in root_chs if c in id_to_idx and y_vals[id_to_idx[c]] is not None]
    if root_ys and root_id in id_to_idx:
        y_vals[id_to_idx[root_id]] = sum(root_ys) / len(root_ys)  # type: ignore[arg-type]
    out["y"] = y_vals
    return out


# Branch midpoints & Angles
def calculate_branch_mid(res: DataFrame, layout: str = "rectangular") -> DataFrame:
    """Exact port of ggtree:::calculate_branch_mid."""
    out = res.copy()
    nodes = [int(n) for n in out["node"].tolist()]
    parents = [int(p) for p in out["parent"].tolist()]
    xs: List[Optional[float]] = [None if _is_scalar_na(v) else float(v) for v in out["x"].tolist()]
    ys: List[Optional[float]] = [None if _is_scalar_na(v) else float(v) for v in out["y"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes)}

    unrooted_layouts = {"equal_angle", "daylight", "ape", "tree_and_leaf", "unrooted"}
    if layout in unrooted_layouts:
        branch_y: List[Optional[float]] = []
        for i, p in enumerate(parents):
            py = ys[id_to_idx[p]] if p in id_to_idx else ys[i]
            yi = ys[i]
            if py is None or yi is None:
                branch_y.append(None)
            else:
                branch_y.append((py + yi) / 2.0)
        out["branch.y"] = branch_y

    branch_x: List[Optional[float]] = []
    for i, p in enumerate(parents):
        px = xs[id_to_idx[p]] if p in id_to_idx else xs[i]
        xi = xs[i]
        if px is None or xi is None:
            branch_x.append(None)
        else:
            branch_x.append((px + xi) / 2.0)
    out["branch"] = branch_x

    if "branch.length" in out:
        br_len = [0.0 if _is_scalar_na(v) else float(v) for v in out["branch.length"].tolist()]
        out["branch.length"] = br_len

    if layout in unrooted_layouts:
        out["branch.x"] = list(branch_x)
    return out


def calculate_angle(data: DataFrame) -> DataFrame:
    """Exact port of ggtree:::calculate_angle."""
    if "y" not in data:
        raise ValueError("Input data must contain 'y' column")
    out = data.copy()
    ys = [float(v) for v in out["y"].tolist() if not _is_scalar_na(v)]
    if not ys:
        out["angle"] = [0.0] * len(out)
        return out
    y_range = max(ys) - min(ys)
    if y_range == 0:
        warnings.warn("Y range is 0, setting all angles to 0")
        out["angle"] = [0.0] * len(out)
    else:
        scale = 360.0 / (y_range + 1.0)
        out["angle"] = [0.0 if _is_scalar_na(v) else scale * float(v) for v in out["y"].tolist()]
    return out


def add_angle_slanted(res: DataFrame) -> DataFrame:
    """Exact port of ggtree:::add_angle_slanted."""
    out = res.copy()
    nodes = [int(n) for n in out["node"].tolist()]
    parents = [int(p) for n, p in zip(nodes, out["parent"].tolist())]
    xs = [float(v) if not _is_scalar_na(v) else 0.0 for v in out["x"].tolist()]
    ys = [float(v) if not _is_scalar_na(v) else 0.0 for v in out["y"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes)}

    rx = (max(xs) - min(xs)) if xs else 1.0
    ry = (max(ys) - min(ys)) if ys else 1.0
    if rx == 0:
        rx = 1.0
    if ry == 0:
        ry = 1.0

    angles = []
    branch_y = []
    for i, p in enumerate(parents):
        pidx = id_to_idx.get(p, i)
        dy = (ys[i] - ys[pidx]) / ry
        dx = (xs[i] - xs[pidx]) / rx
        if dx == 0:
            theta = 0.0
        else:
            theta = math.atan(dy / dx) / math.pi * 180.0
        angles.append(theta)
        branch_y.append((ys[pidx] + ys[i]) / 2.0)

    out["angle"] = angles
    out["branch.y"] = branch_y
    return out


# Unrooted layouts: equal_angle, daylight, ape, tree_and_leaf
def set_branch_length_cladogram(tree: Any) -> Any:
    """Exact port of ggtree:::set_branch_length_cladogram."""
    phylo = as_phylo(tree).copy()
    edge = phylo.edge
    xpos = getXcoord_no_length(phylo)
    phylo.edge_length = [xpos[int(edge[i, 1]) - 1] - xpos[int(edge[i, 0]) - 1] for i in range(edge.shape[0])]
    if isinstance(tree, TreeData):
        td = tree.copy()
        td.phylo = phylo
        return td
    return phylo


def getNodes_by_postorder(tree: Phylo) -> List[int]:
    """Exact port of ggtree:::getNodes_by_postorder."""
    reordered = reorder_phylo(tree, "postorder")
    edge = reordered.edge
    # R: unique(rev(as.vector(t(tree$edge[, c(2, 1)]))))
    flat: List[int] = []
    for i in range(edge.shape[0]):
        flat.append(int(edge[i, 1]))
        flat.append(int(edge[i, 0]))
    seen = set()
    res = []
    for nd in reversed(flat):
        if nd not in seen:
            seen.add(nd)
            res.append(nd)
    return res


def getNodesBreadthFirst_df(df: DataFrame) -> List[int]:
    """Exact port of ggtree:::getNodesBreadthFirst.df."""
    root = getRoot(df)
    nodes = [int(n) for n in df["node"].tolist()]
    parents = [int(p) for p in df["parent"].tolist()]
    kids: Dict[int, List[int]] = {}
    for p, n in zip(parents, nodes):
        if p != n:
            kids.setdefault(p, []).append(n)
    res = [root]
    i = 0
    while i < len(res) and len(res) < len(df):
        cur = res[i]
        i += 1
        if cur in kids:
            res.extend(kids[cur])
    return res


def layoutEqualAngle(model: Any, branch_length: str = "branch.length") -> TblTree:
    """Exact port of ggtree:::layoutEqualAngle."""
    tree = as_phylo(model).copy()
    if tree.edge_length is not None and any(_is_scalar_na(v) for v in tree.edge_length):
        tree.edge_length = None
    if tree.edge_length is None or branch_length == "none":
        tree = set_branch_length_cladogram(tree)

    N = Nnode2(tree)
    brlen = [0.0] * N
    edge = tree.edge
    if tree.edge_length is not None:
        for i in range(edge.shape[0]):
            brlen[int(edge[i, 1]) - 1] = float(tree.edge_length[i])

    root = getRoot(tree)
    df = as_tibble(model)
    parents_set = {int(p) for p, n in zip(df["parent"].tolist(), df["node"].tolist()) if int(p) != int(n)}
    df["isTip"] = [int(n) not in parents_set for n in df["node"].tolist()]
    df = df.sort_values("node")

    xs = [0.0] * N
    ys = [0.0] * N
    starts = [0.0] * N
    ends = [0.0] * N
    angles = [0.0] * N
    starts[root - 1] = 0.0
    ends[root - 1] = 2.0

    # nb.sp[i] = number of tip descendants (self_include = TRUE)
    nb_sp = [0] * N
    ntip = Ntip(tree)
    for i in range(1, ntip + 1):
        nb_sp[i - 1] = 1
    post_edges = reorder_phylo(tree, "postorder").edge
    for i in range(post_edges.shape[0]):
        p, ch = int(post_edges[i, 0]), int(post_edges[i, 1])
        nb_sp[p - 1] += nb_sp[ch - 1]

    children_of: Dict[int, List[int]] = {}
    for i in range(edge.shape[0]):
        p, ch = int(edge[i, 0]), int(edge[i, 1])
        if p != ch:
            children_of.setdefault(p, []).append(ch)

    nodes_order = getNodes_by_postorder(tree)
    for cur_node in nodes_order:
        cur_ntip = nb_sp[cur_node - 1]
        if cur_ntip == 0:
            continue
        chs = children_of.get(cur_node, [])
        st = starts[cur_node - 1]
        en = ends[cur_node - 1]
        cur_x = xs[cur_node - 1]
        cur_y = ys[cur_node - 1]
        total_angle = en - st
        for ch in chs:
            ntip_ch = nb_sp[ch - 1]
            alpha = total_angle * float(ntip_ch) / float(cur_ntip)
            beta = st + alpha / 2.0
            len_ch = brlen[ch - 1]
            xs[ch - 1] = cur_x + math.cos(math.pi * beta) * len_ch
            ys[ch - 1] = cur_y + math.sin(math.pi * beta) * len_ch
            sgn = 1.0 if (beta - 1.0) > 0 else (-1.0 if (beta - 1.0) < 0 else 0.0)
            angles[ch - 1] = -90.0 - 180.0 * beta * sgn
            starts[ch - 1] = st
            ends[ch - 1] = st + alpha
            st += alpha

    df["branch.length"] = [brlen[int(n) - 1] for n in df["node"].tolist()]
    df["x"] = [xs[int(n) - 1] for n in df["node"].tolist()]
    df["y"] = [ys[int(n) - 1] for n in df["node"].tolist()]
    df["start"] = [starts[int(n) - 1] for n in df["node"].tolist()]
    df["end"] = [ends[int(n) - 1] for n in df["node"].tolist()]
    df["angle"] = [angles[int(n) - 1] for n in df["node"].tolist()]
    return TblTree(df, attrs=df.attrs)


def _get_node_angle_pi(x_orig: float, y_orig: float, x: float, y: float) -> float:
    return math.atan2(y - y_orig, x - x_orig) / math.pi


def _get_subtree_unrooted_df(df: DataFrame, node: int) -> List[Dict[str, Any]]:
    """Exact port of ggtree:::getSubtreeUnrooted.df."""
    chs = [int(n) for n in child(df, node)["node"].tolist()]
    if not chs:
        return []
    subtrees: List[Dict[str, Any]] = []
    used_nodes: set = set()
    for c_id in chs:
        st_nodes = [int(n) for n in offspring(df, c_id, self_include=True)["node"].tolist()]
        subtrees.append({"node": c_id, "subtree": st_nodes})
        used_nodes.update(st_nodes)
    all_nodes = [int(n) for n in df["node"].tolist()]
    remaining = [n for n in all_nodes if n not in used_nodes]
    p_df = parent(df, node)
    if len(p_df) > 0 and remaining:
        pid = int(p_df["node"][0])
        subtrees.append({"node": pid, "subtree": remaining})
    return subtrees


def _get_tree_arc_angles(df: DataFrame, origin_id: int, subtree: Dict[str, Any]) -> Tuple[float, float]:
    """Exact port of ggtree:::getTreeArcAngles."""
    df_sorted = df.sort_values("node")
    df_x = [float(v) for v in df_sorted["x"].tolist()]
    df_y = [float(v) for v in df_sorted["y"].tolist()]
    x_origin = df_x[origin_id - 1]
    y_origin = df_y[origin_id - 1]

    subtree_root_id = int(subtree["node"])
    subtree_node_ids = [int(n) for n in subtree["subtree"]]
    origin_children = [int(n) for n in child(df_sorted, origin_id)["node"].tolist()]

    if subtree_root_id in origin_children:
        theta_left = _get_node_angle_pi(x_origin, y_origin, df_x[subtree_root_id - 1], df_y[subtree_root_id - 1])
        theta_right = theta_left
    elif subtree_root_id == origin_id:
        chs = origin_children
        if len(chs) == 2:
            theta1 = _get_node_angle_pi(x_origin, y_origin, df_x[chs[0] - 1], df_y[chs[0] - 1])
            theta2 = _get_node_angle_pi(x_origin, y_origin, df_x[chs[1] - 1], df_y[chs[1] - 1])
            delta = theta1 - theta2
            delta_adj = delta - 2.0 if delta > 1.0 else (delta + 2.0 if delta < -1.0 else delta)
            if delta_adj >= 0:
                theta_left, theta_right = theta1, theta2
            else:
                theta_left, theta_right = theta2, theta1
        else:
            theta_left = _get_node_angle_pi(x_origin, y_origin, df_x[chs[0] - 1], df_y[chs[0] - 1])
            theta_right = theta_left
    else:
        tree_root = getRoot(df_sorted)
        theta_left = _get_node_angle_pi(x_origin, y_origin, df_x[tree_root - 1], df_y[tree_root - 1])
        theta_right = theta_left

    if not subtree_node_ids:
        return (0.0, 0.0)

    arc_left = theta_left
    arc_right = theta_right
    parents_set = set(int(p) for p in df_sorted["parent"].tolist())
    sub_parents = [n for n in subtree_node_ids if n in parents_set and n != origin_id]

    for pid in sub_parents:
        theta_parent = _get_node_angle_pi(x_origin, y_origin, df_x[pid - 1], df_y[pid - 1])
        chs = [int(n) for n in child(df_sorted, pid)["node"].tolist() if int(n) != origin_id]
        for cid in chs:
            theta_child = _get_node_angle_pi(x_origin, y_origin, df_x[cid - 1], df_y[cid - 1])
            if (arc_left < arc_right and not (arc_left < theta_child < arc_right)) or (
                arc_left > arc_right and (arc_right < theta_child < arc_left)
            ):
                pass
            else:
                continue
            delta = theta_child - theta_parent
            delta_adj = delta - 2.0 if delta > 1.0 else (delta + 2.0 if delta < -1.0 else delta)
            theta_child_adj = theta_child
            if delta_adj > 0:
                if abs(delta) > 1.0:
                    if arc_left > 0 and theta_child < 0:
                        theta_child_adj = theta_child + 2.0
                    elif arc_left < 0 and theta_child > 0:
                        theta_child_adj = theta_child - 2.0
                if arc_left < theta_child_adj:
                    arc_left = theta_child
            elif delta_adj < 0:
                if abs(delta) > 1.0:
                    if arc_right > 0 and theta_child < 0:
                        theta_child_adj = theta_child + 2.0
                    elif arc_right < 0 and theta_child > 0:
                        theta_child_adj = theta_child - 2.0
                if arc_right > theta_child_adj:
                    arc_right = theta_child

    if arc_left < 0:
        arc_left += 2.0
    if arc_right < 0:
        arc_right += 2.0
    return (arc_left, arc_right)


def rotateTreePoints_df(df: DataFrame, pivot_node: int, nodes: Sequence[int], angle: float) -> DataFrame:
    """Exact port of ggtree:::rotateTreePoints.df."""
    out = df.sort_values("node")
    cospitheta = math.cos(math.pi * angle)
    sinpitheta = math.sin(math.pi * angle)
    xs = [float(v) for v in out["x"].tolist()]
    ys = [float(v) for v in out["y"].tolist()]
    angles = [float(v) for v in out["angle"].tolist()]
    parents = [int(p) for p in out["parent"].tolist()]

    pivot_x = xs[pivot_node - 1]
    pivot_y = ys[pivot_node - 1]
    node_set = set(int(n) for n in nodes)

    for i in range(len(xs)):
        nid = i + 1
        if nid in node_set:
            dx = xs[i] - pivot_x
            dy = ys[i] - pivot_y
            xs[i] = cospitheta * dx - sinpitheta * dy + pivot_x
            ys[i] = sinpitheta * dx + cospitheta * dy + pivot_y

    parent_set = set(parents)
    tip_rot_set = {n for n in node_set if n not in parent_set}
    for i in range(len(xs)):
        nid = i + 1
        if nid in tip_rot_set:
            p = parents[i]
            ang = _get_node_angle_pi(xs[p - 1], ys[p - 1], xs[i], ys[i])
            angles[i] = 180.0 * (2.0 + ang if ang < 0 else ang)

    out["x"] = xs
    out["y"] = ys
    out["angle"] = angles
    return out


def applyLayoutDaylight(df: DataFrame, node_id: int) -> Tuple[DataFrame, float]:
    """Exact port of ggtree:::applyLayoutDaylight."""
    subtrees = _get_subtree_unrooted_df(df, node_id)
    if len(subtrees) <= 2:
        return (df, 0.0)
    angle_records = []
    for idx, st in enumerate(subtrees):
        left, right = _get_tree_arc_angles(df, node_id, st)
        beta = left - right
        if beta < 0:
            beta += 2.0
        angle_records.append({"left": left, "beta": beta, "subtree_id": idx})
    angle_records.sort(key=lambda r: r["left"])

    total_daylight = 2.0 - sum(r["beta"] for r in angle_records)
    d = total_daylight / len(subtrees)
    new_left_angle = angle_records[0]["left"]
    max_change = 0.0
    out_df = df
    for i in range(1, len(angle_records)):
        new_left_angle = new_left_angle + d + angle_records[i]["beta"]
        adjust_angle = new_left_angle - angle_records[i]["left"]
        max_change = max(max_change, abs(adjust_angle))
        st_id = angle_records[i]["subtree_id"]
        st_nodes = subtrees[st_id]["subtree"]
        out_df = rotateTreePoints_df(out_df, node_id, st_nodes, adjust_angle)
    return (out_df, max_change)


def layoutDaylight(model: Any, branch_length: str = "branch.length", MAX_COUNT: int = 5) -> TblTree:
    """Exact port of ggtree:::layoutDaylight."""
    MINIMUM_AVERAGE_ANGLE_CHANGE = 0.05
    tree_df = layoutEqualAngle(model, branch_length=branch_length)
    bfs_nodes = getNodesBreadthFirst_df(tree_df)
    is_tip_map = {int(n): bool(t) for n, t in zip(tree_df["node"].tolist(), tree_df["isTip"].tolist())}
    internal_bfs = [n for n in bfs_nodes if not is_tip_map.get(n, False)]
    if not internal_bfs:
        return tree_df

    for _ in range(int(MAX_COUNT)):
        total_max = 0.0
        for cur_id in internal_bfs:
            tree_df, mchange = applyLayoutDaylight(tree_df, cur_id)
            total_max += mchange
        ave_change = total_max / len(internal_bfs)
        if ave_change <= MINIMUM_AVERAGE_ANGLE_CHANGE:
            break
    return TblTree(tree_df, attrs=tree_df.attrs)


def layoutApe(model: Any, branch_length: str = "branch.length") -> TblTree:
    """Exact port of ggtree:::layoutApe (using ape::unrooted.xy algorithm)."""
    tree = reorder_phylo(as_phylo(model), "postorder")
    if tree.edge_length is not None and any(_is_scalar_na(v) for v in tree.edge_length):
        tree.edge_length = None
    if tree.edge_length is None or branch_length == "none":
        tree = set_branch_length_cladogram(tree)

    ntip = Ntip(tree)
    nnode = tree.Nnode
    N = ntip + nnode
    edge = tree.edge
    edge_len = tree.edge_length if tree.edge_length is not None else [1.0] * edge.shape[0]

    # ape::node.depth(tree)
    nb_sp = [0.0] * N
    for i in range(ntip):
        nb_sp[i] = 1.0
    for i in range(edge.shape[0]):
        p, ch = int(edge[i, 0]), int(edge[i, 1])
        nb_sp[p - 1] += nb_sp[ch - 1]

    xx = [0.0] * N
    yy = [0.0] * N
    angle = [0.0] * N
    axis = [0.0] * N

    def foo(node: int, ANGLE: float, AXIS: float) -> None:
        inds = [i for i in range(edge.shape[0]) if int(edge[i, 0]) == node]
        sons = [int(edge[i, 1]) for i in inds]
        start = AXIS - ANGLE / 2.0
        for idx_i, son in zip(inds, sons):
            h = float(edge_len[idx_i])
            alpha = ANGLE * nb_sp[son - 1] / nb_sp[node - 1]
            angle[son - 1] = alpha
            beta = start + alpha / 2.0
            axis[son - 1] = beta
            start += alpha
            xx[son - 1] = h * math.cos(beta) + xx[node - 1]
            yy[son - 1] = h * math.sin(beta) + yy[node - 1]
        for son in sons:
            if son > ntip:
                foo(son, angle[son - 1], axis[son - 1])

    foo(ntip + 1, 2.0 * math.pi, 0.0)
    min_xx = min(xx) if xx else 0.0
    min_yy = min(yy) if yy else 0.0

    df = as_tibble(model).sort_values("node")
    parents_set = {int(p) for p, n in zip(df["parent"].tolist(), df["node"].tolist()) if int(p) != int(n)}
    df["isTip"] = [int(n) not in parents_set for n in df["node"].tolist()]
    df["x"] = [xx[int(n) - 1] - min_xx for n in df["node"].tolist()]
    df["y"] = [yy[int(n) - 1] - min_yy for n in df["node"].tolist()]
    return TblTree(df, attrs=df.attrs)


def layoutTreeAndLeaf(model: Any, branch_length: str = "branch.length", niter: int = 150) -> TblTree:
    """
    Port of tree_and_leaf layout: starts from equal_angle / daylight unrooted layout
    and applies force-directed leaf spacing relaxation while preserving tree topology.
    """
    df = layoutEqualAngle(model, branch_length=branch_length).sort_values("node")
    xs = np.array([float(v) for v in df["x"].tolist()], dtype=float)
    ys = np.array([float(v) for v in df["y"].tolist()], dtype=float)
    parents = [int(p) - 1 for p in df["parent"].tolist()]
    is_tip = np.array([bool(t) for t in df["isTip"].tolist()], dtype=bool)
    n = len(xs)

    # Normalize initial coordinates and run spring + leaf-repulsion relaxation
    span = max(float(np.ptp(xs)), float(np.ptp(ys)), 1e-6)
    target_d = span / max(math.sqrt(n), 2.0) * 0.35
    edge_rest = np.array([math.hypot(xs[i] - xs[parents[i]], ys[i] - ys[parents[i]]) for i in range(n)])

    for it in range(int(niter)):
        temp = 0.08 * (1.0 - it / max(1, int(niter)))
        fx = np.zeros(n, dtype=float)
        fy = np.zeros(n, dtype=float)

        # Pairwise short-range repulsion (stronger between tips)
        dx_mat = xs[:, None] - xs[None, :]
        dy_mat = ys[:, None] - ys[None, :]
        dist_mat = np.hypot(dx_mat, dy_mat) + 1e-9
        mask = (dist_mat < target_d * 1.5) & (~np.eye(n, dtype=bool))
        if np.any(mask):
            rep_mag = np.where(mask, (target_d * 1.5 - dist_mat) / dist_mat * 0.15, 0.0)
            fx += np.sum(dx_mat * rep_mag, axis=1)
            fy += np.sum(dy_mat * rep_mag, axis=1)

        # Edge spring force keeping branch lengths close to rest length
        for i in range(n):
            p = parents[i]
            if p == i:
                continue
            dx = xs[i] - xs[p]
            dy = ys[i] - ys[p]
            d = math.hypot(dx, dy) + 1e-9
            diff = (d - edge_rest[i]) / d * 0.4
            fx[i] -= dx * diff
            fy[i] -= dy * diff
            fx[p] += dx * diff * 0.5
            fy[p] += dy * diff * 0.5

        xs += fx * temp
        ys += fy * temp

    df["x"] = xs.tolist()
    df["y"] = ys.tolist()
    return TblTree(df, attrs=df.attrs)


# Non-layered Tidy Tree Layout (exact iterative port of van der Ploeg 2014)
class _VdpNode:
    __slots__ = (
        "id",
        "w",
        "y",
        "c",
        "x",
        "prelim",
        "mod",
        "shift",
        "change",
        "tl",
        "tr",
        "el",
        "er",
        "msel",
        "mser",
    )

    def __init__(self, node_id: int, y_depth: float, width: float = 1.0):
        self.id = node_id
        self.w = width
        self.y = y_depth  # In van der Ploeg, y is depth (= our x coordinate)
        self.c: List["_VdpNode"] = []
        self.x = 0.0
        self.prelim = 0.0
        self.mod = 0.0
        self.shift = 0.0
        self.change = 0.0
        self.tl: Optional["_VdpNode"] = None
        self.tr: Optional["_VdpNode"] = None
        self.el: "_VdpNode" = self
        self.er: "_VdpNode" = self
        self.msel = 0.0
        self.mser = 0.0


class _IYL:
    __slots__ = ("lowY", "index", "nxt")

    def __init__(self, lowY: float, index: int, nxt: Optional["_IYL"]):
        self.lowY = lowY
        self.index = index
        self.nxt = nxt


def _vdp_update_iyl(minY: float, i: int, ih: Optional[_IYL]) -> _IYL:
    while ih is not None and minY >= ih.lowY:
        ih = ih.nxt
    return _IYL(minY, i, ih)


def _vdp_next_left(t: _VdpNode) -> Optional[_VdpNode]:
    return t.tl if not t.c else t.c[0]


def _vdp_next_right(t: _VdpNode) -> Optional[_VdpNode]:
    return t.tr if not t.c else t.c[-1]


def _vdp_set_left_thread(t: _VdpNode, i: int, cl: _VdpNode, modsumcl: float) -> None:
    li = t.c[0].el
    li.tl = cl
    diff = (modsumcl - cl.mod) - t.c[0].msel
    li.mod += diff
    li.prelim -= diff
    t.c[0].el = t.c[i].el
    t.c[0].msel = t.c[i].msel


def _vdp_set_right_thread(t: _VdpNode, i: int, sr: _VdpNode, modsumsr: float) -> None:
    ri = t.c[i].er
    ri.tr = sr
    diff = (modsumsr - sr.mod) - t.c[i].mser
    ri.mod += diff
    ri.prelim -= diff
    t.c[i].er = t.c[i - 1].er
    t.c[i].mser = t.c[i - 1].mser


def _vdp_position_root(t: _VdpNode) -> None:
    t.prelim = (
        t.c[0].prelim
        + t.c[0].mod
        - t.c[0].w / 2.0
        + t.c[-1].prelim
        + t.c[-1].mod
        + t.c[-1].w / 2.0
    ) / 2.0


def _vdp_set_extremes(t: _VdpNode) -> None:
    if not t.c:
        t.el = t
        t.er = t
        t.msel = 0.0
        t.mser = 0.0
    else:
        t.el = t.c[0].el
        t.msel = t.c[0].msel
        t.er = t.c[-1].er
        t.mser = t.c[-1].mser


def _vdp_move_subtree(t: _VdpNode, i: int, si: int, dist: float) -> None:
    t.c[i].mod += dist
    t.c[i].msel += dist
    t.c[i].mser += dist
    if si != i - 1:
        nr = float(i - si)
        t.c[si + 1].shift += dist / nr
        t.c[i].shift -= dist / nr
        t.c[i].change -= dist - dist / nr


def _vdp_add_child_spacing(t: _VdpNode) -> None:
    d = 0.0
    modsumdelta = 0.0
    for ch in t.c:
        d += ch.shift
        modsumdelta += d + ch.change
        ch.mod += modsumdelta


def _vdp_separate(t: _VdpNode, i: int, ih: _IYL) -> None:
    sr: Optional[_VdpNode] = t.c[i - 1]
    mssr = sr.mod
    cl: Optional[_VdpNode] = t.c[i]
    mscl = cl.mod
    cur_ih: Optional[_IYL] = ih
    while sr is not None and cl is not None:
        while cur_ih is not None and cur_ih.nxt is not None and sr.y > cur_ih.lowY:
            cur_ih = cur_ih.nxt
        dist = (mssr + sr.prelim + sr.w / 2.0) - (mscl + cl.prelim - cl.w / 2.0)
        if dist > 0:
            mscl += dist
            _vdp_move_subtree(t, i, cur_ih.index if cur_ih is not None else i - 1, dist)
        sy = sr.y
        cy = cl.y
        if abs(sy - cy) <= 1e-12:
            adv_sr = bool(sr.c) or (not cl.c)
            adv_cl = bool(cl.c) or (not sr.c)
        else:
            adv_sr = sy < cy
            adv_cl = cy < sy
        if adv_sr:
            sr = _vdp_next_right(sr)
            if sr is not None:
                mssr += sr.mod
        if adv_cl:
            cl = _vdp_next_left(cl)
            if cl is not None:
                mscl += cl.mod
    if sr is None and cl is not None:
        _vdp_set_left_thread(t, i, cl, mscl)
    elif sr is not None and cl is None:
        _vdp_set_right_thread(t, i, sr, mssr)


def getYcoord_tidy(tr: Phylo, xpos: Sequence[float]) -> List[float]:
    """
    Iterative non-layered tidy tree layout (van der Ploeg 2014).
    Avoids Python recursion so deep trees (2000+ levels) never overflow the call stack.
    """
    N = getNodeNum(tr)
    root_id = getRoot(tr)
    nodes = {i: _VdpNode(i, float(xpos[i - 1]), width=1.0) for i in range(1, N + 1)}

    edge = tr.edge
    for idx in range(edge.shape[0]):
        p, ch = int(edge[idx, 0]), int(edge[idx, 1])
        if p != ch:
            nodes[p].c.append(nodes[ch])

    # Iterative post-order traversal for firstWalk
    post_order: List[_VdpNode] = []
    stack1 = [nodes[root_id]]
    while stack1:
        u = stack1.pop()
        post_order.append(u)
        for ch in u.c:
            stack1.append(ch)
    post_order.reverse()

    for t in post_order:
        if not t.c:
            _vdp_set_extremes(t)
        else:
            ih = _vdp_update_iyl(t.c[0].el.y, 0, None)
            for i in range(1, len(t.c)):
                minY = t.c[i].er.y
                _vdp_separate(t, i, ih)
                ih = _vdp_update_iyl(minY, i, ih)
            _vdp_position_root(t)
            _vdp_set_extremes(t)

    # Iterative pre-order traversal for secondWalk
    stack2: List[Tuple[_VdpNode, float]] = [(nodes[root_id], 0.0)]
    while stack2:
        t, modsum = stack2.pop()
        modsum += t.mod
        t.x = t.prelim + modsum
        _vdp_add_child_spacing(t)
        for ch in reversed(t.c):
            stack2.append((ch, modsum))

    raw_y = [nodes[i].x for i in range(1, N + 1)]
    min_y = min(raw_y) if raw_y else 0.0
    return [v - min_y for v in raw_y]


# Time & Trait scaling helpers
def scaleX_by_time_from_mrsd(df: DataFrame, mrsd: Any, as_Date: bool = False) -> DataFrame:
    """Exact port of ggtree:::scaleX_by_time_from_mrsd."""
    out = df.copy()
    date_dec = float(Date2decimal(mrsd))
    xs = [float(v) for v in out["x"].tolist()]
    max_x = max(xs) if xs else 0.0
    new_x = [v + date_dec - max_x for v in xs]
    out["x"] = new_x

    nodes = [int(n) for n in out["node"].tolist()]
    parents = [int(p) for p in out["parent"].tolist()]
    id_to_idx = {n: i for i, n in enumerate(nodes)}
    new_branch = [(new_x[id_to_idx.get(p, i)] + new_x[i]) / 2.0 for i, p in enumerate(parents)]
    out["branch"] = new_branch

    if as_Date:
        out["x"] = decimal2Date(new_x)
        out["branch"] = decimal2Date(new_branch)
    return out


def _assign_parent_status(tr: Phylo, df: DataFrame, variable: str) -> DataFrame:
    out = df.sort_values("node")
    yy = list(out[variable].tolist())
    if any(_is_scalar_na(v) for v in yy):
        nodes_post = getNodes_by_postorder(tr)
        for cur in nodes_post:
            chs = child(tr, cur)
            for ch in chs:
                if _is_scalar_na(yy[ch - 1]):
                    yy[ch - 1] = yy[cur - 1]
    out[variable] = yy
    return out


def _assign_child_status(
    tr: Phylo,
    df: DataFrame,
    variable: str,
    yscale_mapping: Optional[Dict[Any, float]] = None,
) -> DataFrame:
    out = df.sort_values("node")
    yy_raw = list(out[variable].tolist())
    if yscale_mapping is not None:
        yy: List[Optional[float]] = [
            float(yscale_mapping[v]) if not _is_scalar_na(v) and v in yscale_mapping else None
            for v in yy_raw
        ]
    else:
        yy = [float(v) if not _is_scalar_na(v) else None for v in yy_raw]

    if any(_is_scalar_na(v) for v in yy):
        nodes_pre = list(reversed(getNodes_by_postorder(tr)))
        for cur in nodes_pre:
            p = parent(tr, cur)
            if p == 0:
                continue
            if _is_scalar_na(yy[p - 1]):
                chs = child(tr, p)
                ch_vals = [float(yy[c - 1]) for c in chs if not _is_scalar_na(yy[c - 1])]
                if ch_vals:
                    yy[p - 1] = sum(ch_vals) / len(ch_vals)
    out[variable] = yy
    return out


def scaleY(
    phylo: Phylo,
    df: DataFrame,
    yscale: str = "none",
    layout: str = "rectangular",
    yscale_mapping: Optional[Dict[Any, float]] = None,
    **kwargs: Any,
) -> DataFrame:
    """Exact port of ggtree:::scaleY."""
    if yscale is None or yscale == "none":
        return df
    if yscale not in df:
        warnings.warn(f"yscale variable '{yscale}' is not available in data frame")
        return df

    out = df.copy()
    vals = [v for v in out[yscale].tolist() if not _is_scalar_na(v)]
    is_num = all(isinstance(v, (int, float, np.number)) and not isinstance(v, bool) for v in vals)
    if is_num:
        tmp = _assign_parent_status(phylo, out, yscale)
        tmp = _assign_child_status(phylo, tmp, yscale)
        y_col = [0.0 if _is_scalar_na(v) else float(v) for v in tmp[yscale].tolist()]
    else:
        if yscale_mapping is None:
            raise ValueError("yscale is categorical; please provide yscale_mapping dict")
        tmp = _assign_parent_status(phylo, out, yscale)
        tmp = _assign_child_status(phylo, tmp, yscale, yscale_mapping=yscale_mapping)
        y_col = [0.0 if _is_scalar_na(v) else float(v) for v in tmp[yscale].tolist()]

    # Align back by node order of out
    node_to_y = {int(n): yv for n, yv in zip(tmp["node"].tolist(), y_col)}
    out["y"] = [node_to_y[int(n)] for n in out["node"].tolist()]
    return out


def adjust_hclust_tip_edge_len(df: DataFrame, phylo: Any, layout: str, branch_length: str = "branch.length") -> DataFrame:
    """Exact port of ggtree:::adjust_hclust_tip_edge_len."""
    p = as_phylo(phylo)
    tip_edge_len = p.attrs.get("tip.edge.len")
    if tip_edge_len is not None and branch_length != "none" and layout in ("rectangular", "dendrogram", "roundrect"):
        out = df.copy()
        xs = [float(v) for v in out["x"].tolist()]
        mx = max(xs) if xs else 0.0
        new_x = [v - mx for v in xs]
        new_br = [float(v) - mx for v in out["branch"].tolist()]
        is_tip = [bool(t) for t in out["isTip"].tolist()]
        nodes = [int(n) for n in out["node"].tolist()]
        for i, (tflag, nid) in enumerate(zip(is_tip, nodes)):
            if tflag and 1 <= nid <= len(tip_edge_len):
                new_x[i] = float(tip_edge_len[nid - 1])
        out["x"] = new_x
        out["branch"] = new_br
        out.attrs["revts.done"] = True
        return out
    return df


# Main tree_layout and fortify functions
def tree_layout(
    model: Any,
    layout: Union[str, Callable[[Any], DataFrame]] = "rectangular",
    ladderize: bool = True,
    right: bool = False,
    branch_length: str = "branch.length",
    mrsd: Any = None,
    as_Date: bool = False,
    yscale: str = "none",
    root_position: float = 0.0,
    **kwargs: Any,
) -> TblTree:
    """
    Compute fortified tree DataFrame with x, y, branch, angle, isTip coordinates
    (exact port of treelayout::tree_layout and ggtree::fortify.phylo).
    """
    branch_length = kwargs.pop("branch.length", branch_length)
    as_Date = kwargs.pop("as.Date", as_Date)
    root_position = kwargs.pop("root.position", root_position)

    if isinstance(model, (MultiPhylo, TreeDataList, list, tuple)):
        return fortify(
            model,
            layout=layout,
            ladderize=ladderize,
            right=right,
            branch_length=branch_length,
            mrsd=mrsd,
            as_Date=as_Date,
            yscale=yscale,
            root_position=root_position,
            **kwargs,
        )

    # Handle callable external layout function
    if callable(layout):
        raw_out = layout(model)
        if isinstance(raw_out, dict) and "xy" in raw_out and not ("x" in raw_out and "y" in raw_out):
            raw_out = raw_out["xy"]
        if is_dataframe_like(raw_out):
            custom_df = ensure_dataframe(raw_out).copy()
        else:
            arr = np.asarray(raw_out)
            if arr.ndim == 2 and arr.shape[1] >= 2:
                custom_df = DataFrame({
                    "node": list(range(1, arr.shape[0] + 1)),
                    "x": arr[:, 0].tolist(),
                    "y": arr[:, 1].tolist(),
                })
            else:
                custom_df = DataFrame(raw_out)
        base_df = as_tibble(model)
        if "node" not in custom_df.columns:
            custom_df["node"] = list(range(1, len(custom_df) + 1))
        parents_set = {int(p) for p, n in zip(base_df["parent"].tolist(), base_df["node"].tolist()) if int(p) != int(n)}
        base_df["isTip"] = [int(n) not in parents_set for n in base_df["node"].tolist()]
        xy_cols = [c for c in ("node", "x", "y") if c in custom_df.columns]
        res = full_join(base_df, custom_df._select_cols(xy_cols), by="node")
        res = calculate_branch_mid(res, layout="external")
        res = calculate_angle(res)
        res.attrs["layout"] = "external"
        return TblTree(res, attrs=res.attrs)

    layout_str = str(layout)
    if layout_str == "unrooted":
        layout_str = "daylight"

    x_phylo = as_phylo(model).copy()
    if ladderize:
        x_phylo = ladderize_phylo(x_phylo, right=right)

    if x_phylo.edge_length is not None and any(_is_scalar_na(v) for v in x_phylo.edge_length):
        warnings.warn("'edge.length' contains NA values; setting 'edge.length' to None automatically.")
        x_phylo.edge_length = None

    # Propagate custom branch_length column from TreeData if requested
    if isinstance(model, TreeData) and branch_length not in ("branch.length", "none"):
        df_td = as_tibble(model)
        if branch_length in df_td:
            val_map = {int(n): v for n, v in zip(df_td["node"].tolist(), df_td[branch_length].tolist())}
            new_el = []
            has_na = False
            for i in range(x_phylo.edge.shape[0]):
                ch = int(x_phylo.edge[i, 1])
                v = val_map.get(ch, None)
                if _is_scalar_na(v):
                    has_na = True
                else:
                    try:
                        new_el.append(float(v))
                    except Exception:
                        has_na = True
            if not has_na and len(new_el) == x_phylo.edge.shape[0]:
                x_phylo.edge_length = new_el
            else:
                x_phylo.edge_length = None

    if layout_str in ("equal_angle", "daylight", "ape", "tree_and_leaf"):
        # Preserve ladderized phylo while keeping model annotations if TreeData
        if isinstance(model, TreeData):
            td_mod = model.copy()
            td_mod.phylo = x_phylo
            target_model: Any = td_mod
        else:
            target_model = x_phylo

        if layout_str == "equal_angle":
            res = layoutEqualAngle(target_model, branch_length=branch_length)
        elif layout_str == "daylight":
            max_cnt = kwargs.get("MAX_COUNT", 5)
            res = layoutDaylight(target_model, branch_length=branch_length, MAX_COUNT=max_cnt)
        elif layout_str == "ape":
            res = layoutApe(target_model, branch_length=branch_length)
        else:
            res = layoutTreeAndLeaf(target_model, branch_length=branch_length)
    else:
        N = Nnode(x_phylo, internal_only=False)
        if x_phylo.edge_length is None or branch_length == "none":
            if layout_str == "slanted":
                sbp, all_nodes = _convert_tips2ancestors_sbp(x_phylo, include_root=True)
                xpos = getXcoord_no_length_slanted(sbp, all_nodes)
                ypos = getYcoord_no_length_slanted(sbp, all_nodes)
            else:
                xpos = getXcoord_no_length(x_phylo)
                ypos = getYcoord_tidy(x_phylo, xpos) if layout_str == "tidy" else getYcoord(x_phylo)
        else:
            xpos = getXcoord(x_phylo)
            ypos = getYcoord_tidy(x_phylo, xpos) if layout_str == "tidy" else getYcoord(x_phylo)

        xypos = DataFrame({
            "node": list(range(1, N + 1)),
            "x": [v + float(root_position) for v in xpos],
            "y": ypos,
        })
        df = as_tibble(model)
        parents_set = {int(p) for p, n in zip(df["parent"].tolist(), df["node"].tolist()) if int(p) != int(n)}
        df["isTip"] = [int(n) not in parents_set for n in df["node"].tolist()]
        res = full_join(df, xypos, by="node")

    res = calculate_branch_mid(res, layout=layout_str)
    if mrsd is not None:
        res = scaleX_by_time_from_mrsd(res, mrsd, as_Date=as_Date)

    if layout_str == "slanted":
        res = add_angle_slanted(res)
    elif layout_str not in ("equal_angle", "daylight"):
        res = calculate_angle(res)

    res = scaleY(x_phylo, res, yscale=yscale, layout=layout_str, **kwargs)
    res = adjust_hclust_tip_edge_len(res, x_phylo, layout_str, branch_length=branch_length)
    res.attrs["layout"] = layout_str
    return TblTree(res, attrs=res.attrs)


def fortify(
    model: Any,
    data: Any = None,
    layout: Union[str, Callable[[Any], DataFrame]] = "rectangular",
    ladderize: bool = True,
    right: bool = False,
    branch_length: str = "branch.length",
    mrsd: Any = None,
    as_Date: bool = False,
    yscale: str = "none",
    root_position: float = 0.0,
    **kwargs: Any,
) -> TblTree:
    """
    Fortify a phylogenetic tree (Phylo, TreeData, MultiPhylo, TreeDataList, Tanglegram,
    Hclust, or DataFrame) into a TblTree DataFrame ready for plotting
    (exact port of R/method-fortify.R).
    """
    branch_length = kwargs.pop("branch.length", branch_length)
    as_Date = kwargs.pop("as.Date", as_Date)
    root_position = kwargs.pop("root.position", root_position)

    if is_dataframe_like(model):
        model_df = ensure_dataframe(model)
        if "x" in model_df and "y" in model_df:
            out = TblTree(model_df.copy(), attrs=model_df.attrs)
            if "layout" not in out.attrs:
                out.attrs["layout"] = layout if isinstance(layout, str) else "external"
            return out
        return tree_layout(
            as_phylo(model_df),
            layout=layout,
            ladderize=ladderize,
            right=right,
            branch_length=branch_length,
            mrsd=mrsd,
            as_Date=as_Date,
            yscale=yscale,
            root_position=root_position,
            **kwargs,
        )

    # Tanglegram (gtangle) support
    if hasattr(model, "is_gtangle") and getattr(model, "is_gtangle", False):
        return model.fortify(
            layout=layout,
            ladderize=ladderize,
            right=right,
            branch_length=branch_length,
            mrsd=mrsd,
            as_Date=as_Date,
            yscale=yscale,
            **kwargs,
        )

    # MultiPhylo / TreeDataList / list of trees
    if isinstance(model, (MultiPhylo, TreeDataList, list, tuple)):
        dfs = [
            fortify(
                tr,
                layout=layout,
                ladderize=ladderize,
                right=right,
                branch_length=branch_length,
                mrsd=mrsd,
                as_Date=as_Date,
                yscale=yscale,
                root_position=root_position,
                **kwargs,
            )
            for tr in model
        ]
        names = getattr(model, "names", None)
        if names is None:
            names = [f"Tree #{i+1}" for i in range(len(dfs))]
        for d, nm in zip(dfs, names):
            d[".id"] = nm
        combined = bind_rows(*dfs)
        from .dataframe import factor
        combined[".id"] = factor(combined[".id"].tolist(), levels=list(names))
        combined.attrs["layout"] = layout if isinstance(layout, str) else "external"
        return TblTree(combined, attrs=combined.attrs)

    if isinstance(model, Hclust):
        hang = kwargs.pop("hang", 0.1)
        phy = model.to_phylo(hang=hang)
        return tree_layout(
            phy,
            layout=layout,
            ladderize=ladderize,
            right=right,
            branch_length=branch_length,
            mrsd=mrsd,
            as_Date=as_Date,
            yscale=yscale,
            root_position=root_position,
            **kwargs,
        )

    return tree_layout(
        model,
        layout=layout,
        ladderize=ladderize,
        right=right,
        branch_length=branch_length,
        mrsd=mrsd,
        as_Date=as_Date,
        yscale=yscale,
        root_position=root_position,
        **kwargs,
    )


# Layout attribute & tip alignment helpers (from R/layout.R)
def get_layout(plot: Any) -> Optional[str]:
    """Retrieve layout string from a ggtree plot or fortified DataFrame (returns None if unknown)."""
    if plot is None:
        return None
    if is_dataframe_like(plot):
        return ensure_dataframe(plot).attrs.get("layout", None)
    if hasattr(plot, "data") and is_dataframe_like(plot.data) and "layout" in ensure_dataframe(plot.data).attrs:
        return ensure_dataframe(plot.data).attrs["layout"]
    if hasattr(plot, "layout") and isinstance(plot.layout, str) and plot.layout:
        return plot.layout
    return None


def tip_y_alignment(x: Any, tol: float = 1e-6) -> str:
    """
    Inspect tip y coordinates and return 'ok', 'duplicated', or 'uneven'
    (exact port of tip_y_alignment in R/utilities.R).
    """
    raw_df = x.data if hasattr(x, "data") and is_dataframe_like(x.data) else x
    if not is_dataframe_like(raw_df):
        return "ok"
    df = ensure_dataframe(raw_df)
    if "y" not in df or "isTip" not in df:
        return "ok"
    tips_df = df[df["isTip"]]
    y = sorted(float(v) for v in tips_df["y"].tolist() if not _is_scalar_na(v))
    if len(y) <= 1:
        return "ok"
    dy = [y[i + 1] - y[i] for i in range(len(y) - 1)]
    if any(d <= tol for d in dy):
        return "duplicated"
    if (max(dy) - min(dy)) > tol:
        return "uneven"
    return "ok"


def adjust_cladelabel_angle(angle: float, horizontal: bool = True) -> float:
    """Port of ggtree:::adjust_cladelabel_angle in R/geom_cladelabel.R."""
    ang = float(angle)
    if horizontal:
        return ang
    if ang > 180.0:
        return ang + 90.0
    return ang + 270.0


# Clade & Strip Geometry Helpers (used by geom_hilight, geom_cladelab, geom_strip, geom_balance)
def check_reverse(data: DataFrame) -> bool:
    """Check whether tree x coordinates are reversed (tips have smaller x than internal nodes)."""
    tips = data[data["isTip"]]
    internals = data[~data["isTip"]]
    if len(tips) == 0 or len(internals) == 0:
        return False
    mean_tip_x = tips["x"].mean()
    mean_int_x = internals["x"].mean()
    return bool(mean_tip_x < mean_int_x)


def get_clade_position_(data: DataFrame, node: Union[int, str], reverse: bool = False) -> DataFrame:
    """Exact port of ggtree:::get_clade_position_."""
    nid = taxa2node(data, node)
    sp_df = offspring(data, nid, self_include=True)
    if len(sp_df) == 0:
        sp_df = data[data["node"] == nid]
    xs = [float(v) for v in sp_df["x"].tolist() if not _is_scalar_na(v)]
    ys = [float(v) for v in sp_df["y"].tolist() if not _is_scalar_na(v)]
    node_row = data[data["node"] == nid]
    br_len = 0.0
    if "branch.length" in data and len(node_row) > 0 and not _is_scalar_na(node_row["branch.length"][0]):
        br_len = float(node_row["branch.length"][0])

    if "branch.length" in data:
        if reverse:
            xmin = min(xs)
            xmax = max(xs) + br_len / 2.0
        else:
            xmin = min(xs) - br_len / 2.0
            xmax = max(xs)
    else:
        brs = [float(v) for v in sp_df["branch"].tolist() if not _is_scalar_na(v)] if "branch" in sp_df else xs
        xmin = min(brs)
        xmax = max(xs)

    w = float(getOption("clade_width_extend", 0.5))
    return DataFrame({
        "xmin": [xmin],
        "xmax": [xmax],
        "ymin": [min(ys) - w],
        "ymax": [max(ys) + w],
    })


def get_clade_position(treeview: Any, node: Union[int, str]) -> DataFrame:
    """
    Get bounding box (xmin, xmax, ymin, ymax) of clade rooted at 'node'
    (exact port of ggtree::get_clade_position in R/geom_hilight.R).
    """
    df = treeview.data if hasattr(treeview, "data") else treeview
    rev = check_reverse(df)
    return get_clade_position_(df, node, reverse=rev)


def clade_geometry(
    trdf: DataFrame,
    nodeids: Sequence[Union[int, str]],
    offset: Union[float, Sequence[float]] = 0.0,
    align: Union[bool, Sequence[bool]] = False,
    extend: Union[float, Sequence[float]] = 0.0,
    angle: Union[str, float, Sequence[Any]] = "auto",
    horizontal: Union[bool, Sequence[bool]] = True,
    mode: str = "rectangular",
    label: Optional[Sequence[Any]] = None,
) -> Dict[str, DataFrame]:
    """
    Compute bar and label position DataFrames for geom_cladelab / geom_cladelabel
    (exact port of clade_geometry / build_cladebar_df / build_cladelabel_df).
    """
    n = len(nodeids)
    offsets = list(offset) if isinstance(offset, (list, tuple, np.ndarray, ColumnSeries)) else [float(offset)] * n
    aligns = list(align) if isinstance(align, (list, tuple, np.ndarray, ColumnSeries)) else [bool(align)] * n
    if isinstance(extend, (list, tuple, np.ndarray)) and len(extend) == 2 and n != 2:
        extends = [tuple(extend)] * n
    elif isinstance(extend, (list, tuple, np.ndarray, ColumnSeries)):
        extends = [(float(e), float(e)) if not isinstance(e, (list, tuple)) else (float(e[0]), float(e[1])) for e in extend]
    else:
        extends = [(float(extend), float(extend))] * n
    angles = list(angle) if isinstance(angle, (list, tuple, np.ndarray, ColumnSeries)) and not isinstance(angle, str) else [angle] * n
    horizontals = list(horizontal) if isinstance(horizontal, (list, tuple, np.ndarray, ColumnSeries)) else [bool(horizontal)] * n
    labels = list(label) if label is not None else [""] * n

    max_tree_x = max((float(v) for v in trdf["x"].tolist() if not _is_scalar_na(v)), default=0.0)
    bar_rows = []
    label_rows = []

    for nid_raw, off, aln, ext, ang, horiz, lbl in zip(nodeids, offsets, aligns, extends, angles, horizontals, labels):
        nid = taxa2node(trdf, nid_raw)
        sp_df = offspring(trdf, nid, self_include=True)
        if len(sp_df) == 0:
            sp_df = trdf[trdf["node"] == nid]
        ys = [float(v) for v in sp_df["y"].tolist() if not _is_scalar_na(v)]
        xs = [float(v) for v in sp_df["x"].tolist() if not _is_scalar_na(v)]
        sp_angles = [float(v) for v in sp_df["angle"].tolist() if not _is_scalar_na(v)] if "angle" in sp_df else [0.0]

        mx = max_tree_x if aln else (max(xs) if xs else 0.0)
        y_min = (min(ys) if ys else 0.0) - ext[1]
        y_max = (max(ys) if ys else 0.0) + ext[0]

        if ang == "auto":
            raw_ang = (min(sp_angles) + max(sp_angles)) / 2.0 if sp_angles else 0.0
            calc_ang = adjust_cladelabel_angle(raw_ang, horizontal=bool(horiz))
        else:
            calc_ang = float(ang)

        bx = mx * 1.02 + float(off)
        lx = mx * 1.03 + float(off)
        bar_rows.append({"x": bx, "xend": bx, "y": y_min, "yend": y_max, "angle": 0.0, "node": nid})
        label_rows.append({
            "x": lx,
            "xend": lx,
            "y": (y_min + y_max) / 2.0,
            "yend": y_max,
            "angle": calc_ang,
            "label": lbl,
            "node": nid,
        })

    return {"bar": DataFrame(bar_rows), "label": DataFrame(label_rows)}
