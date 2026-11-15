"""
R-compatible DataFrame, TblTree, ColumnSeries, Factor, and vector/table utilities
for ggtree in Python. Works out-of-the-box with only numpy + Python stdlib,
and interoperates seamlessly with pandas.DataFrame when installed.
"""

from __future__ import annotations

import copy
import csv
import math
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import numpy as np


# Global R-style options store (e.g., options(ignore.negative.edge=True), getOption(...))
_OPTIONS: Dict[str, Any] = {
    "ignore.negative.edge": False,
    "clade_width_extend": 0.5,
    "clade_align": False,
}


def options(**kwargs: Any) -> Dict[str, Any]:
    """Set or query R-style global options. Supports both '_' and '.' in option names."""
    old = dict(_OPTIONS)
    for k, v in kwargs.items():
        _OPTIONS[k] = v
        _OPTIONS[k.replace("_", ".")] = v
        _OPTIONS[k.replace(".", "_")] = v
    return old


def getOption(name: str, default: Any = None) -> Any:
    """Retrieve an R-style global option."""
    if name in _OPTIONS:
        return _OPTIONS[name]
    alt = name.replace("_", ".")
    if alt in _OPTIONS:
        return _OPTIONS[alt]
    return default


def attr(obj: Any, which: str) -> Any:
    """Get R-style attribute from an object."""
    if hasattr(obj, "attrs") and isinstance(obj.attrs, dict):
        if which in obj.attrs:
            return obj.attrs[which]
        alt = which.replace("_", ".")
        if alt in obj.attrs:
            return obj.attrs[alt]
        alt2 = which.replace(".", "_")
        if alt2 in obj.attrs:
            return obj.attrs[alt2]
    if hasattr(obj, f"_{which}"):
        return getattr(obj, f"_{which}")
    return getattr(obj, which, None)


def set_attr(obj: Any, which: str, value: Any) -> Any:
    """Set R-style attribute on an object."""
    if hasattr(obj, "attrs") and isinstance(obj.attrs, dict):
        obj.attrs[which] = value
    else:
        setattr(obj, which, value)
    return obj


# R constants
NA = None
TRUE = True
FALSE = False
NULL = None
LETTERS = [chr(ord("A") + i) for i in range(26)]
letters = [chr(ord("a") + i) for i in range(26)]

# Shared random generator for reproducible set.seed / set_seed
_RNG = np.random.default_rng(42)


def set_seed(seed: Union[int, float]) -> None:
    """Set random seed (R's set.seed)."""
    global _RNG
    s = int(abs(round(float(seed))))
    np.random.seed(s)
    _RNG = np.random.default_rng(s)


def _get_rng() -> np.random.Generator:
    """Return the active shared random generator updated by set_seed()."""
    return _RNG


def is_na(val: Any) -> Union[bool, "ColumnSeries", np.ndarray]:
    """Check whether value(s) are NA / None / NaN / pd.NA / pd.NaT."""
    if isinstance(val, ColumnSeries):
        return val.isna()
    if _is_series_like(val):
        return ColumnSeries(val).isna()
    if isinstance(val, (list, tuple, np.ndarray, Factor)):
        return ColumnSeries(list(val)).isna()
    if val is None or type(val).__name__ in ("NAType", "NaTType"):
        return True
    if isinstance(val, np.datetime64):
        return bool(np.isnat(val))
    if isinstance(val, (float, np.floating)):
        return bool(math.isnan(float(val)))
    try:
        neq = val != val
        if type(neq).__name__ == "NAType":
            return True
        if neq:  # NaN check
            return True
    except Exception:
        pass
    return False


def _is_scalar_na(val: Any) -> bool:
    if val is None or type(val).__name__ in ("NAType", "NaTType"):
        return True
    if isinstance(val, np.datetime64):
        return bool(np.isnat(val))
    if isinstance(val, (float, np.floating)):
        return bool(math.isnan(float(val)))
    if isinstance(val, str) and val == "<NA>":
        return True
    try:
        neq = val != val
        if type(neq).__name__ == "NAType":
            return True
        return bool(neq)
    except Exception:
        return False


def _clean_external_scalar(val: Any) -> Any:
    """Normalize pandas/polars/pyarrow scalar NA sentinels and numpy scalars."""
    if val is None or type(val).__name__ in ("NAType", "NaTType"):
        return None
    if isinstance(val, np.datetime64):
        return None if np.isnat(val) else val.item()
    if isinstance(val, np.generic):
        if isinstance(val, np.floating) and math.isnan(float(val)):
            return float("nan")
        return val.item()
    return val


def is_dataframe_like(obj: Any) -> bool:
    """
    Return True for ggtree.DataFrame or external DataFrame objects
    (pandas, polars, pyarrow.Table/RecordBatch, modin, cudf, geopandas, or __dataframe__ interchange).
    Does not require external libraries to be installed.
    """
    if isinstance(obj, DataFrame):
        return True
    if obj is None or isinstance(obj, (str, bytes, dict, list, tuple, np.ndarray, ColumnSeries, Factor)):
        return False
    if hasattr(obj, "layers") or type(obj).__name__ in (
        "Phylo",
        "TreeData",
        "MultiPhylo",
        "TreeDataList",
        "Hclust",
    ):
        return False
    # pandas / modin / cudf / geopandas DataFrame
    if hasattr(obj, "columns") and hasattr(obj, "to_dict"):
        return True
    # polars DataFrame
    if hasattr(obj, "columns") and hasattr(obj, "to_dicts") and hasattr(obj, "schema"):
        return True
    # pyarrow Table / RecordBatch
    if hasattr(obj, "column_names") and hasattr(obj, "to_pydict"):
        return True
    # Python DataFrame Interchange Protocol (__dataframe__)
    df_proto = getattr(obj, "__dataframe__", None)
    if callable(df_proto):
        return True
    return False


def _is_series_like(obj: Any) -> bool:
    """Return True for external 1D Series / Index / Categorical / Array objects."""
    if obj is None or isinstance(obj, (DataFrame, ColumnSeries, Factor, str, bytes, dict, list, tuple, np.ndarray)):
        return False
    if is_dataframe_like(obj):
        return False
    if hasattr(obj, "layers") or type(obj).__name__ in (
        "Phylo",
        "TreeData",
        "MultiPhylo",
        "TreeDataList",
        "Hclust",
    ):
        return False
    tname = type(obj).__name__
    if tname in (
        "Series",
        "Index",
        "RangeIndex",
        "CategoricalIndex",
        "DatetimeIndex",
        "MultiIndex",
        "Int64Index",
        "Float64Index",
        "StringIndex",
        "Categorical",
        "Array",
        "ChunkedArray",
        "IntegerArray",
        "FloatingArray",
        "BooleanArray",
        "StringArray",
        "ArrowExtensionArray",
    ):
        return True
    if hasattr(obj, "to_list") and hasattr(obj, "dtype") and not hasattr(obj, "columns"):
        return True
    if hasattr(obj, "to_pylist") and not hasattr(obj, "column_names"):
        return True
    if hasattr(obj, "tolist") and hasattr(obj, "dtype") and getattr(obj, "ndim", 1) == 1 and not hasattr(obj, "columns"):
        return True
    return False


def _series_to_list_and_levels(s: Any) -> Tuple[List[Any], Optional[List[Any]]]:
    """Extract Python value list and optional ordered category levels from a Series/Categorical/Index."""
    levels: Optional[List[Any]] = None
    try:
        if hasattr(s, "categories") and hasattr(s, "codes"):
            levels = [_clean_external_scalar(x) for x in list(s.categories)]
        else:
            cat = getattr(s, "cat", None)
            if cat is not None:
                if hasattr(cat, "categories"):
                    levels = [_clean_external_scalar(x) for x in list(cat.categories)]
                elif hasattr(cat, "get_categories"):
                    cats_obj = cat.get_categories()
                    raw_cats = cats_obj.to_list() if hasattr(cats_obj, "to_list") else list(cats_obj)
                    levels = [_clean_external_scalar(x) for x in raw_cats]
    except Exception:
        levels = None

    if hasattr(s, "tolist"):
        raw = s.tolist()
    elif hasattr(s, "to_list"):
        raw = s.to_list()
    elif hasattr(s, "to_pylist"):
        raw = s.to_pylist()
    else:
        raw = list(s)

    vals = [
        None if (levels is not None and _is_scalar_na(v)) else _clean_external_scalar(v)
        for v in raw
    ]
    return vals, levels


def ensure_dataframe(obj: Any) -> "DataFrame":
    """Return obj if already a ggtree.DataFrame, else convert external DataFrame-like input."""
    if isinstance(obj, DataFrame):
        return obj
    if obj is None:
        return DataFrame()
    return DataFrame(obj)


class Factor:
    """R-style categorical factor with ordered levels."""

    def __init__(
        self,
        values: Iterable[Any],
        levels: Optional[Sequence[Any]] = None,
        labels: Optional[Sequence[Any]] = None,
        ordered: bool = False,
    ):
        ext_levels: Optional[List[Any]] = None
        if _is_series_like(values):
            raw, ext_levels = _series_to_list_and_levels(values)
            if hasattr(values, "cat") and getattr(values.cat, "ordered", False):
                ordered = ordered or bool(values.cat.ordered)
            elif hasattr(values, "ordered") and bool(getattr(values, "ordered", False)):
                ordered = ordered or bool(values.ordered)
        else:
            raw = [_clean_external_scalar(v) for v in values]

        if levels is None:
            if ext_levels is not None:
                self.levels = list(ext_levels)
            else:
                seen = []
                for v in raw:
                    if not _is_scalar_na(v) and v not in seen:
                        seen.append(v)
                # Sort levels naturally as strings/numbers like R's factor()
                try:
                    seen = sorted(seen, key=lambda x: (0, float(x)) if str(x).replace(".", "", 1).replace("-", "", 1).isdigit() else (1, str(x)))
                except Exception:
                    seen = sorted(seen, key=str)
                self.levels = [str(x) if not isinstance(x, (int, float, str)) else x for x in seen]
        else:
            self.levels = list(levels)

        if labels is not None:
            label_map = dict(zip(self.levels, labels))
            self.values = [None if _is_scalar_na(v) else label_map.get(v, None) for v in raw]
            self.levels = list(labels)
        else:
            lvl_set = set(self.levels) | {str(x) for x in self.levels}
            self.values = [
                None if _is_scalar_na(v) or (v not in lvl_set and str(v) not in lvl_set) else v
                for v in raw
            ]
        self.ordered = ordered

    def __len__(self) -> int:
        return len(self.values)

    def __iter__(self) -> Iterator[Any]:
        return iter(self.values)

    def __getitem__(self, idx: Any) -> Any:
        if isinstance(idx, (int, np.integer)):
            return self.values[int(idx)]
        s = ColumnSeries(self.values, levels=self.levels)[idx]
        return Factor(list(s), levels=self.levels, ordered=self.ordered)

    def __eq__(self, other: Any) -> "ColumnSeries":  # type: ignore[override]
        return ColumnSeries(self.values, levels=self.levels) == other

    def __ne__(self, other: Any) -> "ColumnSeries":  # type: ignore[override]
        return ColumnSeries(self.values, levels=self.levels) != other

    def tolist(self) -> List[Any]:
        return list(self.values)

    def __repr__(self) -> str:
        return f"Factor({self.values!r}, levels={self.levels!r})"


def factor(
    x: Any,
    levels: Optional[Sequence[Any]] = None,
    labels: Optional[Sequence[Any]] = None,
    ordered: bool = False,
) -> Factor:
    """Create an R-style Factor."""
    if isinstance(x, (int, float, str, bool)):
        x = [x]
    return Factor(x, levels=levels, labels=labels, ordered=ordered)


class ColumnSeries:
    """1D column vector supporting vectorized arithmetic, comparisons, and boolean indexing."""

    def __init__(
        self,
        data: Iterable[Any],
        name: Optional[str] = None,
        levels: Optional[List[Any]] = None,
    ):
        if isinstance(data, ColumnSeries):
            self._data = list(data._data)
            self.name = name if name is not None else data.name
            self.levels = list(levels) if levels is not None else (list(data.levels) if data.levels is not None else None)
        elif isinstance(data, Factor):
            self._data = list(data.values)
            self.name = name
            self.levels = list(levels) if levels is not None else list(data.levels)
        elif _is_series_like(data):
            vals, cat_lvls = _series_to_list_and_levels(data)
            self._data = vals
            raw_name = getattr(data, "name", None)
            self.name = name if name is not None else (str(raw_name) if raw_name is not None else None)
            self.levels = list(levels) if levels is not None else cat_lvls
        elif isinstance(data, (str, bytes)):
            self._data = [data]
            self.name = name
            self.levels = list(levels) if levels is not None else None
        elif isinstance(data, Iterable):
            self._data = [_clean_external_scalar(v) for v in data]
            self.name = name
            self.levels = list(levels) if levels is not None else None
        else:
            self._data = [_clean_external_scalar(data)]
            self.name = name
            self.levels = list(levels) if levels is not None else None

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self) -> Iterator[Any]:
        return iter(self._data)

    def tolist(self) -> List[Any]:
        return list(self._data)

    def to_numpy(self) -> np.ndarray:
        try:
            return np.array(
                [np.nan if _is_scalar_na(v) else v for v in self._data],
                dtype=float,
            )
        except Exception:
            return np.array(self._data, dtype=object)

    @property
    def values(self) -> np.ndarray:
        return self.to_numpy()

    def copy(self) -> "ColumnSeries":
        return ColumnSeries(list(self._data), name=self.name, levels=self.levels)

    def __getitem__(self, idx: Any) -> Any:
        if isinstance(idx, (int, np.integer)):
            return self._data[int(idx)]
        if isinstance(idx, slice):
            return ColumnSeries(self._data[idx], name=self.name, levels=self.levels)
        if isinstance(idx, ColumnSeries):
            idx = idx._data
        if isinstance(idx, (list, tuple, np.ndarray)):
            if len(idx) == 0:
                return ColumnSeries([], name=self.name, levels=self.levels)
            if all(isinstance(x, (bool, np.bool_)) or x is None for x in idx):
                return ColumnSeries(
                    [v for v, flag in zip(self._data, idx) if bool(flag)],
                    name=self.name,
                    levels=self.levels,
                )
            return ColumnSeries(
                [self._data[int(i)] if not _is_scalar_na(i) else None for i in idx],
                name=self.name,
                levels=self.levels,
            )
        raise TypeError(f"Unsupported index type: {type(idx)}")

    def __setitem__(self, idx: Any, value: Any) -> None:
        if isinstance(idx, (int, np.integer)):
            self._data[int(idx)] = value
            return
        if isinstance(idx, ColumnSeries):
            idx = idx._data
        if isinstance(idx, (list, tuple, np.ndarray)):
            if len(idx) > 0 and all(isinstance(x, (bool, np.bool_)) or x is None for x in idx):
                if isinstance(value, (ColumnSeries, list, tuple, np.ndarray)) and not isinstance(value, (str, bytes)):
                    v_iter = iter(value)
                    for i, flag in enumerate(idx):
                        if bool(flag):
                            self._data[i] = next(v_iter)
                else:
                    for i, flag in enumerate(idx):
                        if bool(flag):
                            self._data[i] = value
                return
            if isinstance(value, (ColumnSeries, list, tuple, np.ndarray)) and not isinstance(value, (str, bytes)):
                for i, v in zip(idx, value):
                    self._data[int(i)] = v
            else:
                for i in idx:
                    self._data[int(i)] = value
            return
        if isinstance(idx, slice):
            self._data[idx] = value

    def isna(self) -> "ColumnSeries":
        return ColumnSeries([_is_scalar_na(v) for v in self._data], name=self.name)

    def notna(self) -> "ColumnSeries":
        return ColumnSeries([not _is_scalar_na(v) for v in self._data], name=self.name)

    def isin(self, values: Iterable[Any]) -> "ColumnSeries":
        val_set = set(values)
        str_set = {str(v) for v in val_set if not _is_scalar_na(v)}
        res = []
        for v in self._data:
            if _is_scalar_na(v):
                res.append(False)
            else:
                res.append(v in val_set or str(v) in str_set)
        return ColumnSeries(res, name=self.name)

    def unique(self) -> List[Any]:
        seen = []
        for v in self._data:
            if not _is_scalar_na(v) and v not in seen:
                seen.append(v)
        return seen

    def min(self, na_rm: bool = True) -> Any:
        vals = [v for v in self._data if not _is_scalar_na(v)]
        if not vals:
            return np.nan
        return min(vals)

    def max(self, na_rm: bool = True) -> Any:
        vals = [v for v in self._data if not _is_scalar_na(v)]
        if not vals:
            return np.nan
        return max(vals)

    def mean(self, na_rm: bool = True) -> float:
        vals = [float(v) for v in self._data if not _is_scalar_na(v)]
        if not vals:
            return float("nan")
        return sum(vals) / len(vals)

    def sum(self, na_rm: bool = True) -> float:
        vals = [float(v) for v in self._data if not _is_scalar_na(v)]
        return sum(vals)

    def any(self) -> bool:
        return any(bool(v) for v in self._data if not _is_scalar_na(v))

    def all(self) -> bool:
        vals = [bool(v) for v in self._data if not _is_scalar_na(v)]
        return all(vals) if vals else True

    def astype(self, dtype: Any) -> "ColumnSeries":
        res = []
        for v in self._data:
            if _is_scalar_na(v):
                res.append(None)
            elif dtype in (int, np.int64, "int"):
                res.append(int(float(v)))
            elif dtype in (float, np.float64, "float"):
                res.append(float(v))
            elif dtype in (str, "str"):
                res.append(str(v))
            elif dtype in (bool, "bool"):
                res.append(bool(v))
            else:
                res.append(dtype(v))
        return ColumnSeries(res, name=self.name)

    def _binop(self, other: Any, op: Callable[[Any, Any], Any], rev: bool = False) -> "ColumnSeries":
        if isinstance(other, ColumnSeries):
            other_list = other._data
        elif isinstance(other, Factor):
            other_list = other.values
        elif isinstance(other, (list, tuple, np.ndarray)):
            other_list = list(other)
        else:
            other_list = [other] * len(self._data)

        out = []
        for a, b in zip(self._data, other_list):
            if _is_scalar_na(a) or _is_scalar_na(b):
                out.append(None)
            else:
                try:
                    out.append(op(b, a) if rev else op(a, b))
                except Exception:
                    out.append(None)
        return ColumnSeries(out, name=self.name)

    def _cmpop(self, other: Any, op: Callable[[Any, Any], bool]) -> "ColumnSeries":
        if isinstance(other, ColumnSeries):
            other_list = other._data
        elif isinstance(other, Factor):
            other_list = other.values
        elif isinstance(other, (list, tuple, np.ndarray)):
            other_list = list(other)
        else:
            other_list = [other] * len(self._data)

        out = []
        for a, b in zip(self._data, other_list):
            if _is_scalar_na(a) or _is_scalar_na(b):
                out.append(False)
            else:
                try:
                    if isinstance(a, (int, float, np.number)) and isinstance(b, str):
                        out.append(op(str(a), b) if not b.replace(".", "", 1).replace("-", "", 1).isdigit() else op(float(a), float(b)))
                    elif isinstance(b, (int, float, np.number)) and isinstance(a, str):
                        out.append(op(a, str(b)) if not a.replace(".", "", 1).replace("-", "", 1).isdigit() else op(float(a), float(b)))
                    else:
                        out.append(bool(op(a, b)))
                except Exception:
                    out.append(False)
        return ColumnSeries(out, name=self.name)

    def __add__(self, other: Any) -> "ColumnSeries":
        return self._binop(other, lambda a, b: a + b)

    def __radd__(self, other: Any) -> "ColumnSeries":
        return self._binop(other, lambda a, b: a + b, rev=True)

    def __sub__(self, other: Any) -> "ColumnSeries":
        return self._binop(other, lambda a, b: a - b)

    def __rsub__(self, other: Any) -> "ColumnSeries":
        return self._binop(other, lambda a, b: a - b, rev=True)

    def __mul__(self, other: Any) -> "ColumnSeries":
        return self._binop(other, lambda a, b: a * b)

    def __rmul__(self, other: Any) -> "ColumnSeries":
        return self._binop(other, lambda a, b: a * b, rev=True)

    def __truediv__(self, other: Any) -> "ColumnSeries":
        return self._binop(other, lambda a, b: a / b)

    def __rtruediv__(self, other: Any) -> "ColumnSeries":
        return self._binop(other, lambda a, b: a / b, rev=True)

    def __neg__(self) -> "ColumnSeries":
        return ColumnSeries([None if _is_scalar_na(v) else -v for v in self._data], name=self.name)

    def __abs__(self) -> "ColumnSeries":
        return ColumnSeries([None if _is_scalar_na(v) else abs(v) for v in self._data], name=self.name)

    def __eq__(self, other: Any) -> "ColumnSeries":  # type: ignore[override]
        if self.levels is not None:
            return self._cmpop(other, lambda a, b: (a == b) or (str(a) == str(b)))
        return self._cmpop(other, lambda a, b: a == b)

    def __ne__(self, other: Any) -> "ColumnSeries":  # type: ignore[override]
        if self.levels is not None:
            return self._cmpop(other, lambda a, b: not ((a == b) or (str(a) == str(b))))
        return self._cmpop(other, lambda a, b: a != b)

    def __lt__(self, other: Any) -> "ColumnSeries":
        return self._cmpop(other, lambda a, b: a < b)

    def __le__(self, other: Any) -> "ColumnSeries":
        return self._cmpop(other, lambda a, b: a <= b)

    def __gt__(self, other: Any) -> "ColumnSeries":
        return self._cmpop(other, lambda a, b: a > b)

    def __ge__(self, other: Any) -> "ColumnSeries":
        return self._cmpop(other, lambda a, b: a >= b)

    def __and__(self, other: Any) -> "ColumnSeries":
        if isinstance(other, ColumnSeries):
            other_list = other._data
        elif isinstance(other, (list, tuple, np.ndarray)):
            other_list = list(other)
        else:
            other_list = [other] * len(self._data)
        return ColumnSeries(
            [bool(a) and bool(b) if not (_is_scalar_na(a) or _is_scalar_na(b)) else False for a, b in zip(self._data, other_list)],
            name=self.name,
        )

    def __or__(self, other: Any) -> "ColumnSeries":
        if isinstance(other, ColumnSeries):
            other_list = other._data
        elif isinstance(other, (list, tuple, np.ndarray)):
            other_list = list(other)
        else:
            other_list = [other] * len(self._data)
        return ColumnSeries(
            [bool(a) or bool(b) if not (_is_scalar_na(a) and _is_scalar_na(b)) else False for a, b in zip(self._data, other_list)],
            name=self.name,
        )

    def __invert__(self) -> "ColumnSeries":
        return ColumnSeries([False if _is_scalar_na(v) else not bool(v) for v in self._data], name=self.name)

    def __repr__(self) -> str:
        return f"ColumnSeries({self._data!r}, name={self.name!r})"


class _LocIndexer:
    def __init__(self, df: "DataFrame", is_iloc: bool = False):
        self._df = df
        self._is_iloc = is_iloc

    def __getitem__(self, key: Any) -> Any:
        if isinstance(key, tuple) and len(key) == 2:
            rows, cols = key
            sub = self._df._slice_rows(rows, is_iloc=self._is_iloc)
            if isinstance(cols, str):
                col_name = sub._resolve_col(cols)
                if isinstance(rows, (int, np.integer)):
                    return sub._cols[col_name][0]
                return sub[col_name]
            if isinstance(cols, (int, np.integer)):
                col_name = sub.columns[int(cols)]
                if isinstance(rows, (int, np.integer)):
                    return sub._cols[col_name][0]
                return sub[col_name]
            if isinstance(cols, (list, tuple)):
                resolved = [sub._resolve_col(c) if isinstance(c, str) else sub.columns[int(c)] for c in cols]
                return sub._select_cols(resolved)
            return sub
        return self._df._slice_rows(key, is_iloc=self._is_iloc)

    def __setitem__(self, key: Any, value: Any) -> None:
        if isinstance(key, tuple) and len(key) == 2:
            rows, cols = key
            row_indices = self._df._resolve_row_indices(rows, is_iloc=self._is_iloc)
            if isinstance(cols, str):
                col_name = self._df._resolve_col(cols, create=True)
                if col_name not in self._df._cols:
                    self._df._cols[col_name] = [None] * len(self._df)
                if isinstance(value, (ColumnSeries, list, tuple, np.ndarray)) and not isinstance(value, (str, bytes)):
                    for r, v in zip(row_indices, value):
                        self._df._cols[col_name][r] = v
                else:
                    for r in row_indices:
                        self._df._cols[col_name][r] = value
            elif isinstance(cols, (list, tuple)):
                for c in cols:
                    self.__setitem__((rows, c), value)


class DataFrame:
    """
    Lightweight R tibble / data.frame and pandas.DataFrame compatible table.
    Supports both R-style and pandas-style indexing, dot-column resolution,
    row names, and arbitrary R attributes (.attrs).
    """

    def __init__(
        self,
        data: Any = None,
        columns: Optional[Sequence[str]] = None,
        index: Optional[Sequence[Any]] = None,
        attrs: Optional[Dict[str, Any]] = None,
        rownames: Optional[Sequence[Any]] = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_cols", {})
        object.__setattr__(self, "_levels", {})
        object.__setattr__(self, "_rownames", None)
        object.__setattr__(self, "attrs", dict(attrs) if attrs else {})
        if rownames is None:
            rownames = kwargs.pop("row.names", None)
        if rownames is not None and index is None:
            index = rownames

        if isinstance(data, DataFrame):
            for k, v in data._cols.items():
                self._cols[k] = list(v)
            self._levels = copy.deepcopy(data._levels)
            self._rownames = list(data._rownames) if data._rownames is not None else None
            self.attrs = copy.deepcopy(data.attrs)
        elif data is not None and hasattr(data, "to_dict") and hasattr(data, "columns"):
            # Convert pandas / modin / cudf / geopandas DataFrame seamlessly
            has_iloc = hasattr(data, "iloc")
            for j, col in enumerate(list(data.columns)):
                cname = str(col)
                try:
                    series = data.iloc[:, j] if has_iloc else data[col]
                except Exception:
                    series = data[col]
                vals, lvls = _series_to_list_and_levels(series)
                self._cols[cname] = vals
                if lvls is not None:
                    self._levels[cname] = lvls
            idx_obj = getattr(data, "index", None)
            if idx_obj is not None:
                is_default_range = (
                    type(idx_obj).__name__ == "RangeIndex"
                    and getattr(idx_obj, "start", None) == 0
                    and getattr(idx_obj, "step", None) == 1
                    and getattr(idx_obj, "name", None) is None
                )
                if not is_default_range:
                    raw_idx = idx_obj.tolist() if hasattr(idx_obj, "tolist") else list(idx_obj)
                    self._rownames = [str(_clean_external_scalar(x)) for x in raw_idx]
            if hasattr(data, "attrs") and isinstance(data.attrs, dict):
                self.attrs.update(copy.deepcopy(data.attrs))
        elif data is not None and hasattr(data, "columns") and hasattr(data, "to_dicts") and hasattr(data, "schema"):
            # Convert polars.DataFrame seamlessly
            for col in list(data.columns):
                cname = str(col)
                vals, lvls = _series_to_list_and_levels(data[col])
                self._cols[cname] = vals
                if lvls is not None:
                    self._levels[cname] = lvls
        elif data is not None and hasattr(data, "column_names") and hasattr(data, "to_pydict"):
            # Convert pyarrow.Table / RecordBatch seamlessly
            pyd = data.to_pydict()
            for col in list(data.column_names):
                cname = str(col)
                self._cols[cname] = [_clean_external_scalar(v) for v in pyd.get(col, [])]
        elif data is not None and callable(getattr(data, "__dataframe__", None)):
            # Convert any Python DataFrame Interchange Protocol object
            converted = False
            try:
                import pandas as pd  # type: ignore

                pdf = pd.api.interchange.from_dataframe(data)
                tmp = DataFrame(pdf)
                self._cols = tmp._cols
                self._levels = tmp._levels
                self._rownames = tmp._rownames
                self.attrs.update(tmp.attrs)
                converted = True
            except Exception:
                pass
            if not converted:
                try:
                    import pyarrow.interchange as pai  # type: ignore

                    pa_tbl = pai.from_dataframe(data)
                    tmp = DataFrame(pa_tbl)
                    self._cols = tmp._cols
                    self._levels = tmp._levels
                    converted = True
                except Exception:
                    pass
            if not converted:
                if hasattr(data, "to_dict"):
                    d_raw = data.to_dict()
                    if isinstance(d_raw, dict):
                        self._init_from_dict({str(k): v for k, v in d_raw.items()})
                        converted = True
                if not converted:
                    raise TypeError(
                        f"Cannot convert __dataframe__ object of type {type(data).__name__} without pandas or pyarrow installed."
                    )
        elif isinstance(data, dict):
            self._init_from_dict(data)
        elif isinstance(data, (list, tuple, np.ndarray)):
            if len(data) == 0:
                if columns:
                    for c in columns:
                        self._cols[str(c)] = []
            elif isinstance(data[0], dict):
                keys = list(columns) if columns else []
                if not keys:
                    for row in data:
                        for k in row.keys():
                            if k not in keys:
                                keys.append(k)
                for k in keys:
                    self._cols[str(k)] = [row.get(k, None) for row in data]
            else:
                arr = np.asarray(data)
                if arr.ndim == 1:
                    col_name = columns[0] if columns else "V1"
                    self._cols[str(col_name)] = arr.tolist()
                else:
                    ncol = arr.shape[1]
                    col_names = list(columns) if columns else [f"V{i+1}" for i in range(ncol)]
                    for j, cname in enumerate(col_names):
                        self._cols[str(cname)] = arr[:, j].tolist()

        if kwargs:
            self._init_from_dict(kwargs)

        if index is not None:
            object.__setattr__(self, "_rownames", [str(x) for x in index])

    def _init_from_dict(self, d: Dict[str, Any]) -> None:
        # Determine target length for scalar recycling (just like R's data.frame / tibble)
        max_len = len(self) if self._cols else 0
        processed: Dict[str, List[Any]] = {}
        for k, v in d.items():
            if isinstance(v, Factor):
                self._levels[k] = list(v.levels)
                lst = list(v.values)
            elif isinstance(v, ColumnSeries):
                if v.levels is not None:
                    self._levels[k] = list(v.levels)
                lst = list(v._data)
            elif _is_series_like(v):
                vals, cat_lvls = _series_to_list_and_levels(v)
                if cat_lvls is not None:
                    self._levels[k] = cat_lvls
                lst = vals
            elif isinstance(v, (list, tuple, np.ndarray)):
                raw_lst = v.tolist() if isinstance(v, np.ndarray) else v
                lst = [_clean_external_scalar(item) for item in raw_lst]
            else:
                lst = [_clean_external_scalar(v)]
            processed[k] = lst
            if len(lst) > max_len:
                max_len = len(lst)

        for k, lst in processed.items():
            if len(lst) == 1 and max_len > 1:
                self._cols[k] = lst * max_len
            elif len(lst) == 0 and max_len == 0:
                self._cols[k] = []
            elif len(lst) < max_len and len(lst) > 0:
                # Recycle vector like R's rep_len
                reps = (max_len + len(lst) - 1) // len(lst)
                self._cols[k] = (lst * reps)[:max_len]
            else:
                self._cols[k] = lst

    def _resolve_col(self, name: str, create: bool = False) -> str:
        if name in self._cols:
            return name
        alt_dot = name.replace("_", ".")
        if alt_dot in self._cols:
            return alt_dot
        alt_under = name.replace(".", "_")
        if alt_under in self._cols:
            return alt_under
        if create:
            return name
        raise KeyError(f"Column '{name}' not found in {list(self._cols.keys())}")

    @property
    def columns(self) -> List[str]:
        return list(self._cols.keys())

    @columns.setter
    def columns(self, new_cols: Sequence[str]) -> None:
        new_cols_list = [str(c) for c in new_cols]
        if len(new_cols_list) != len(self._cols):
            raise ValueError("Length mismatch when renaming columns")
        new_dict = {}
        new_levels = {}
        for old_k, new_k in zip(list(self._cols.keys()), new_cols_list):
            new_dict[new_k] = self._cols[old_k]
            if old_k in self._levels:
                new_levels[new_k] = self._levels[old_k]
        object.__setattr__(self, "_cols", new_dict)
        object.__setattr__(self, "_levels", new_levels)

    @property
    def rownames(self) -> List[str]:
        if self._rownames is not None and len(self._rownames) == len(self):
            return list(self._rownames)
        return [str(i + 1) for i in range(len(self))]

    @rownames.setter
    def rownames(self, names: Optional[Sequence[Any]]) -> None:
        if names is None:
            object.__setattr__(self, "_rownames", None)
        else:
            object.__setattr__(self, "_rownames", [str(x) for x in names])

    @property
    def index(self) -> List[str]:
        return self.rownames

    @index.setter
    def index(self, names: Optional[Sequence[Any]]) -> None:
        self.rownames = names

    def __len__(self) -> int:
        if not self._cols:
            return 0
        first_col = next(iter(self._cols.values()))
        return len(first_col)

    @property
    def shape(self) -> Tuple[int, int]:
        return (len(self), len(self._cols))

    @property
    def empty(self) -> bool:
        return len(self) == 0

    @property
    def loc(self) -> _LocIndexer:
        return _LocIndexer(self, is_iloc=False)

    @property
    def iloc(self) -> _LocIndexer:
        return _LocIndexer(self, is_iloc=True)

    def copy(self) -> "DataFrame":
        cls = type(self)
        out = cls()
        for k, v in self._cols.items():
            out._cols[k] = list(v)
        out._levels = copy.deepcopy(self._levels)
        out._rownames = list(self._rownames) if self._rownames is not None else None
        out.attrs = copy.deepcopy(self.attrs)
        return out

    def _select_cols(self, col_names: Sequence[str]) -> "DataFrame":
        cls = type(self)
        out = cls()
        for c in col_names:
            rc = self._resolve_col(c)
            out._cols[rc] = list(self._cols[rc])
            if rc in self._levels:
                out._levels[rc] = list(self._levels[rc])
        out._rownames = list(self._rownames) if self._rownames is not None else None
        out.attrs = copy.deepcopy(self.attrs)
        return out

    def _resolve_row_indices(self, rows: Any, is_iloc: bool = False) -> List[int]:
        n = len(self)
        if isinstance(rows, slice):
            return list(range(*rows.indices(n)))
        if isinstance(rows, (int, np.integer)):
            r = int(rows)
            if r < 0:
                r += n
            return [r]
        if isinstance(rows, ColumnSeries):
            rows = rows._data
        if isinstance(rows, (list, tuple, np.ndarray)):
            if len(rows) == 0:
                return []
            if all(isinstance(x, (bool, np.bool_)) or x is None for x in rows):
                return [i for i, flag in enumerate(rows) if bool(flag)]
            if all(isinstance(x, str) for x in rows) and self._rownames is not None:
                rmap = {rn: i for i, rn in enumerate(self._rownames)}
                return [rmap[x] for x in rows if x in rmap]
            return [int(x) for x in rows if not _is_scalar_na(x)]
        if isinstance(rows, str) and self._rownames is not None:
            return [i for i, rn in enumerate(self._rownames) if rn == rows]
        raise TypeError(f"Unsupported row indexer: {type(rows)}")

    def _slice_rows(self, rows: Any, is_iloc: bool = False) -> "DataFrame":
        indices = self._resolve_row_indices(rows, is_iloc=is_iloc)
        cls = type(self)
        out = cls()
        for k, v in self._cols.items():
            out._cols[k] = [v[i] if 0 <= i < len(v) else None for i in indices]
        out._levels = copy.deepcopy(self._levels)
        if self._rownames is not None:
            out._rownames = [self._rownames[i] if 0 <= i < len(self._rownames) else str(i + 1) for i in indices]
        out.attrs = copy.deepcopy(self.attrs)
        return out

    def __getitem__(self, key: Any) -> Any:
        # 1. Single column name
        if isinstance(key, str):
            col_name = self._resolve_col(key)
            return ColumnSeries(self._cols[col_name], name=col_name, levels=self._levels.get(col_name))
        # 2. Tuple (rows, cols) - R-style df[rows, cols]
        if isinstance(key, tuple) and len(key) == 2:
            rows, cols = key
            sub = self if (isinstance(rows, slice) and rows == slice(None)) else self._slice_rows(rows)
            if isinstance(cols, slice) and cols == slice(None):
                return sub
            if isinstance(cols, str):
                col_name = sub._resolve_col(cols)
                if isinstance(rows, (int, np.integer)):
                    return sub._cols[col_name][0]
                return ColumnSeries(sub._cols[col_name], name=col_name, levels=sub._levels.get(col_name))
            if isinstance(cols, (int, np.integer)):
                col_name = sub.columns[int(cols)]
                if isinstance(rows, (int, np.integer)):
                    return sub._cols[col_name][0]
                return ColumnSeries(sub._cols[col_name], name=col_name, levels=sub._levels.get(col_name))
            if isinstance(cols, (list, tuple)):
                resolved = [sub._resolve_col(c) if isinstance(c, str) else sub.columns[int(c)] for c in cols]
                return sub._select_cols(resolved)
            return sub
        # 3. Boolean mask or row indices or list of column names
        if isinstance(key, ColumnSeries):
            return self._slice_rows(key._data)
        if isinstance(key, (list, tuple, np.ndarray)):
            if len(key) > 0 and all(isinstance(x, str) for x in key):
                return self._select_cols(key)
            return self._slice_rows(key)
        if isinstance(key, (int, np.integer, slice)):
            return self._slice_rows(key)
        raise KeyError(f"Unsupported key: {key!r}")

    def __setitem__(self, key: Any, value: Any) -> None:
        if isinstance(key, tuple) and len(key) == 2:
            self.loc[key] = value
            return
        if isinstance(key, str):
            col_name = self._resolve_col(key, create=True)
            if value is None:
                # In R, df$col <- NULL deletes the column if df already has rows, or sets NAs?
                # Wait: tidytree does `.data[[group_name]] <- NULL` to delete a column!
                if col_name in self._cols:
                    del self._cols[col_name]
                    self._levels.pop(col_name, None)
                return
            n = len(self)
            if isinstance(value, Factor):
                self._levels[col_name] = list(value.levels)
                lst = list(value.values)
            elif isinstance(value, ColumnSeries):
                if value.levels is not None:
                    self._levels[col_name] = list(value.levels)
                lst = list(value._data)
            elif _is_series_like(value):
                vals, cat_lvls = _series_to_list_and_levels(value)
                if cat_lvls is not None:
                    self._levels[col_name] = cat_lvls
                lst = vals
            elif isinstance(value, (list, tuple, np.ndarray)):
                raw_lst = value.tolist() if isinstance(value, np.ndarray) else value
                lst = [_clean_external_scalar(item) for item in raw_lst]
            else:
                lst = [_clean_external_scalar(value)] * max(1, n)
            if n > 0 and len(lst) == 1 and n > 1:
                lst = lst * n
            self._cols[col_name] = lst
            return
        raise TypeError(f"Unsupported item assignment key: {key!r}")

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_") or name in ("attrs", "columns", "shape", "loc", "iloc", "rownames", "index", "empty"):
            raise AttributeError(name)
        cols = object.__getattribute__(self, "_cols")
        if name in cols:
            levels = object.__getattribute__(self, "_levels")
            return ColumnSeries(cols[name], name=name, levels=levels.get(name))
        alt_dot = name.replace("_", ".")
        if alt_dot in cols:
            levels = object.__getattribute__(self, "_levels")
            return ColumnSeries(cols[alt_dot], name=alt_dot, levels=levels.get(alt_dot))
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith("_") or name in ("attrs", "columns", "rownames", "index"):
            object.__setattr__(self, name, value)
            return
        cols = object.__getattribute__(self, "_cols")
        alt_dot = name.replace("_", ".")
        if alt_dot in cols and name not in cols:
            self[alt_dot] = value
        else:
            self[name] = value

    def __contains__(self, item: str) -> bool:
        try:
            self._resolve_col(item)
            return True
        except KeyError:
            return False

    def get(self, key: str, default: Any = None) -> Any:
        try:
            return self[key]
        except KeyError:
            return default

    def iterrows(self) -> Iterator[Tuple[int, Dict[str, Any]]]:
        n = len(self)
        cols = list(self._cols.keys())
        for i in range(n):
            yield i, {c: self._cols[c][i] for c in cols}

    def to_records(self) -> List[Dict[str, Any]]:
        n = len(self)
        cols = list(self._cols.keys())
        return [{c: self._cols[c][i] for c in cols} for i in range(n)]

    def row_dict(self, idx: int) -> Dict[str, Any]:
        return {c: self._cols[c][int(idx)] for c in self._cols.keys()}

    def to_dict(self, orient: str = "list") -> Any:
        if orient == "records":
            return self.to_records()
        return {k: list(v) for k, v in self._cols.items()}

    def to_pandas(self, categorical: bool = False) -> Any:
        """Convert to a pandas.DataFrame (optionally preserving Factor levels as pd.Categorical)."""
        try:
            import pandas as pd  # type: ignore
        except ImportError as exc:
            raise ImportError(
                "pandas is required for DataFrame.to_pandas(). Install with `pip install pandas`."
            ) from exc
        df = pd.DataFrame(self.to_dict())
        if categorical and self._levels:
            for col, lvls in self._levels.items():
                if col in df.columns and lvls:
                    vals = self._cols[col]
                    lvl_set = set(lvls)
                    vals_conv: List[Any] = []
                    for v in vals:
                        if _is_scalar_na(v):
                            vals_conv.append(None)
                        elif v in lvl_set:
                            vals_conv.append(v)
                        elif str(v) in lvl_set:
                            vals_conv.append(str(v))
                        else:
                            vals_conv.append(v)
                    df[col] = pd.Categorical(vals_conv, categories=list(lvls), ordered=False)
        if self._rownames is not None and len(self._rownames) == len(self):
            df.index = self._rownames
        if self.attrs:
            df.attrs.update(self.attrs)
        return df

    def to_polars(self) -> Any:
        """Convert to a polars.DataFrame."""
        try:
            import polars as pl  # type: ignore
        except ImportError as exc:
            raise ImportError(
                "polars is required for DataFrame.to_polars(). Install with `pip install polars`."
            ) from exc
        return pl.DataFrame(self.to_dict(), strict=False)

    def to_arrow(self) -> Any:
        """Convert to a pyarrow.Table."""
        try:
            import pyarrow as pa  # type: ignore
        except ImportError as exc:
            raise ImportError(
                "pyarrow is required for DataFrame.to_arrow(). Install with `pip install pyarrow`."
            ) from exc
        return pa.table(self.to_dict())

    @classmethod
    def from_pandas(cls, df: Any) -> "DataFrame":
        """Create a ggtree DataFrame from a pandas.DataFrame."""
        return cls(df)

    @classmethod
    def from_polars(cls, df: Any) -> "DataFrame":
        """Create a ggtree DataFrame from a polars.DataFrame."""
        return cls(df)

    @classmethod
    def from_arrow(cls, table: Any) -> "DataFrame":
        """Create a ggtree DataFrame from a pyarrow.Table."""
        return cls(table)

    def head(self, n: int = 6) -> "DataFrame":
        return self._slice_rows(slice(0, n))

    def tail(self, n: int = 6) -> "DataFrame":
        return self._slice_rows(slice(max(0, len(self) - n), len(self)))

    def reset_index(self, drop: bool = True) -> "DataFrame":
        out = self.copy()
        if drop:
            out._rownames = None
        return out

    def sort_values(self, by: Union[str, Sequence[str]], ascending: bool = True) -> "DataFrame":
        if isinstance(by, str):
            by = [by]
        by_cols = [self._resolve_col(c) for c in by]
        n = len(self)

        def sort_key(i: int) -> Tuple[Any, ...]:
            key_parts = []
            for c in by_cols:
                v = self._cols[c][i]
                if _is_scalar_na(v):
                    key_parts.append((1, ""))
                elif isinstance(v, (int, float, np.number)):
                    key_parts.append((0, float(v) if ascending else -float(v)))
                else:
                    key_parts.append((0, str(v)))
            return tuple(key_parts)

        indices = sorted(range(n), key=sort_key)
        if not ascending and all(not isinstance(self._cols[by_cols[0]][i], (int, float, np.number)) for i in range(n)):
            indices.reverse()
        return self._slice_rows(indices)

    def mutate(self, **kwargs: Any) -> "DataFrame":
        return mutate(self, **kwargs)

    def filter(self, *args: Any, **kwargs: Any) -> "DataFrame":
        return filter(self, *args, **kwargs)

    def select(self, *cols: str) -> "DataFrame":
        return select(self, *cols)

    def __repr__(self) -> str:
        n = len(self)
        cols = self.columns
        header = f"# A {type(self).__name__}: {n} x {len(cols)}\n"
        if n == 0 or not cols:
            return header + f"Columns: {cols}"
        show_n = min(n, 15)
        widths = {c: max(len(c), max(len(str(self._cols[c][i])) for i in range(show_n))) for c in cols}
        lines = ["  ".join(c.rjust(widths[c]) for c in cols)]
        for i in range(show_n):
            lines.append("  ".join(str(self._cols[c][i]).rjust(widths[c]) for c in cols))
        if n > show_n:
            lines.append(f"... with {n - show_n} more rows")
        return header + "\n".join(lines)


class TblTree(DataFrame):
    """R's tbl_tree subclass of DataFrame representing a fortified or tabular phylogenetic tree."""
    pass


# R-compatible constructor & vector functions
def data_frame(*args: Any, row_names: Optional[Sequence[Any]] = None, **kwargs: Any) -> DataFrame:
    """R-compatible data.frame(...) constructor."""
    row_names = kwargs.pop("row.names", row_names)
    df = DataFrame(*args, **kwargs)
    if row_names is not None:
        df.rownames = row_names
    return df


def tibble(*args: Any, **kwargs: Any) -> TblTree:
    """R-compatible tibble(...) constructor."""
    return TblTree(*args, **kwargs)


def c(*args: Any, **kwargs: Any) -> Any:
    """R-style vector concatenation c(...). Flattens sequences and returns dict for named vectors."""
    if kwargs and not args:
        return dict(kwargs)
    out: List[Any] = []
    for a in args:
        if isinstance(a, (list, tuple, np.ndarray, ColumnSeries, Factor)) or _is_series_like(a):
            vals = _series_to_list_and_levels(a)[0] if _is_series_like(a) else a
            for item in vals:
                out.append(item)
        else:
            out.append(a)
    if kwargs:
        merged = {str(i): v for i, v in enumerate(out)}
        merged.update(kwargs)
        return merged
    return out


def rep(x: Any, times: Optional[Union[int, Sequence[int]]] = None, each: int = 1, length_out: Optional[int] = None) -> List[Any]:
    """R-style rep(x, times=..., each=..., length_out=...)."""
    vec = list(x) if isinstance(x, (list, tuple, np.ndarray, ColumnSeries, Factor)) and not isinstance(x, (str, bytes)) else [x]
    if each > 1:
        vec = [item for item in vec for _ in range(int(each))]
    if isinstance(times, (list, tuple, np.ndarray)):
        out = []
        for item, t in zip(vec, times):
            out.extend([item] * int(t))
        vec = out
    elif times is not None:
        vec = vec * int(times)
    if length_out is not None:
        if len(vec) == 0:
            return []
        reps = (int(length_out) + len(vec) - 1) // len(vec)
        vec = (vec * reps)[: int(length_out)]
    return vec


def seq(from_: float = 1, to: float = 1, by: Optional[float] = None, length_out: Optional[int] = None) -> List[Union[int, float]]:
    """R-style seq(from, to, by=..., length_out=...)."""
    if length_out is not None:
        arr = np.linspace(from_, to, int(length_out))
        return arr.tolist()
    if by is None:
        by = 1 if to >= from_ else -1
    if isinstance(from_, int) and isinstance(to, int) and isinstance(by, int):
        step = int(by)
        end = int(to) + (1 if step > 0 else -1)
        return list(range(int(from_), end, step))
    res = []
    cur = float(from_)
    if by > 0:
        while cur <= to + 1e-12:
            res.append(cur)
            cur += by
    else:
        while cur >= to - 1e-12:
            res.append(cur)
            cur += by
    return res


def seq_len(n: int) -> List[int]:
    """R's seq_len(n): 1..n."""
    return list(range(1, int(n) + 1))


def seq_along(x: Sequence[Any]) -> List[int]:
    """R's seq_along(x): 1..len(x)."""
    return list(range(1, len(x) + 1))


def paste(*args: Any, sep: str = " ", collapse: Optional[str] = None) -> Union[str, List[str]]:
    """R's paste(..., sep=' ', collapse=None)."""
    vecs = [
        list(a) if isinstance(a, (list, tuple, np.ndarray, ColumnSeries, Factor)) and not isinstance(a, (str, bytes)) else [a]
        for a in args
    ]
    if not vecs:
        return "" if collapse is not None else []
    max_len = max(len(v) for v in vecs)
    if max_len == 0:
        return "" if collapse is not None else []
    recycled = [(v * ((max_len + len(v) - 1) // len(v)))[:max_len] for v in vecs]
    joined = [sep.join(str(recycled[j][i]) for j in range(len(vecs))) for i in range(max_len)]
    if collapse is not None:
        return collapse.join(joined)
    if max_len == 1 and all(not isinstance(a, (list, tuple, np.ndarray, ColumnSeries, Factor)) for a in args):
        return joined[0]
    return joined


def paste0(*args: Any, collapse: Optional[str] = None) -> Union[str, List[str]]:
    """R's paste0(..., collapse=None)."""
    return paste(*args, sep="", collapse=collapse)


def rnorm(n: int, mean: float = 0.0, sd: float = 1.0) -> np.ndarray:
    """R's rnorm(n, mean=0, sd=1)."""
    return _RNG.normal(loc=mean, scale=sd, size=int(n))


def runif(n: int, min: float = 0.0, max: float = 1.0) -> np.ndarray:
    """R's runif(n, min=0, max=1)."""
    return _RNG.uniform(low=min, high=max, size=int(n))


def rbinom(n: int, size: int = 1, prob: float = 0.5) -> np.ndarray:
    """R's rbinom(n, size, prob)."""
    return _RNG.binomial(n=int(size), p=float(prob), size=int(n))


def sample(x: Sequence[Any], size: Optional[int] = None, replace: bool = False) -> List[Any]:
    """R's sample(x, size, replace=FALSE)."""
    if isinstance(x, int):
        pool = list(range(1, x + 1))
    else:
        pool = list(x)
    k = len(pool) if size is None else int(size)
    idx = _RNG.choice(len(pool), size=k, replace=replace)
    return [pool[int(i)] for i in idx]


def match(x: Iterable[Any], table: Iterable[Any], nomatch: Any = None) -> List[Any]:
    """0-based index match helper (or None if not found)."""
    tbl = list(table)
    lookup: Dict[Any, int] = {}
    for idx, val in enumerate(tbl):
        if val not in lookup:
            lookup[val] = idx
    return [lookup.get(v, nomatch) for v in x]


# dplyr-style verbs
def mutate(df: Any, _extra_env: Optional[Dict[str, Any]] = None, **kwargs: Any) -> DataFrame:
    """R dplyr::mutate equivalent."""
    import inspect
    from .aes import SymExpr, eval_expr
    extra_env = dict(_extra_env or {})
    if not extra_env:
        frame = inspect.currentframe()
        if frame is not None and frame.f_back is not None:
            extra_env.update(frame.f_back.f_globals)
            extra_env.update(frame.f_back.f_locals)
    out = ensure_dataframe(df).copy()
    for k, expr in kwargs.items():
        if isinstance(expr, SymExpr):
            out[k] = expr.eval(out)
        elif callable(expr):
            out[k] = expr(out)
        elif isinstance(expr, str) and expr in out:
            out[k] = out[expr]
        elif isinstance(expr, str):
            val = eval_expr(out, expr, extra_env=extra_env)
            out[k] = val
        else:
            out[k] = expr
    return out


def filter(df: Any, *conditions: Any, _extra_env: Optional[Dict[str, Any]] = None, **kwargs: Any) -> DataFrame:
    """R dplyr::filter equivalent."""
    import inspect
    from .aes import SymExpr, eval_subset
    extra_env = dict(_extra_env or {})
    if not extra_env:
        frame = inspect.currentframe()
        if frame is not None and frame.f_back is not None:
            extra_env.update(frame.f_back.f_globals)
            extra_env.update(frame.f_back.f_locals)
    out = ensure_dataframe(df).copy()
    for cond in conditions:
        if isinstance(cond, SymExpr):
            mask = eval_subset(out, cond, extra_env=extra_env)
        elif callable(cond):
            mask = cond(out)
        elif isinstance(cond, str):
            mask = eval_subset(out, cond, extra_env=extra_env)
        else:
            mask = cond
        out = out[mask]
    for k, v in kwargs.items():
        out = out[out[k] == v]
    return out


def select(df: Any, *cols: str) -> DataFrame:
    """R dplyr::select equivalent."""
    return ensure_dataframe(df)._select_cols(cols)


def arrange(df: Any, *cols: str, descending: bool = False) -> DataFrame:
    """R dplyr::arrange equivalent."""
    return ensure_dataframe(df).sort_values(list(cols), ascending=not descending)


def left_join(x: Any, y: DataFrame, by: Optional[Union[str, Sequence[str], Dict[str, str]]] = None) -> Any:
    """R dplyr::left_join / tidytree::left_join."""
    from .tree import Phylo, TreeData, as_treedata, as_tibble
    if isinstance(x, (Phylo, TreeData)):
        td = as_treedata(x)
        tbl = as_tibble(td)
        df_y = y if isinstance(y, DataFrame) else DataFrame(y)
        joined = _join_dataframes(tbl, df_y, by=by, how="left")
        core_cols = {"parent", "branch.length", "label"}
        extra_cols = ["node"] + [c for c in joined.columns if c not in core_cols and c != "node"]
        td.data = joined._select_cols(extra_cols) if len(extra_cols) > 1 else DataFrame()
        return td

    df_x = x if isinstance(x, DataFrame) else DataFrame(x)
    df_y = y if isinstance(y, DataFrame) else DataFrame(y)
    return _join_dataframes(df_x, df_y, by=by, how="left")


def full_join(x: Any, y: DataFrame, by: Optional[Union[str, Sequence[str], Dict[str, str]]] = None) -> Any:
    """R dplyr::full_join / tidytree::full_join."""
    from .tree import Phylo, TreeData, as_treedata, as_tibble
    if isinstance(x, (Phylo, TreeData)):
        td = as_treedata(x)
        tbl = as_tibble(td)
        df_y = y if isinstance(y, DataFrame) else DataFrame(y)
        joined = _join_dataframes(tbl, df_y, by=by, how="full")
        if "node" in joined.columns:
            joined = joined[[not _is_scalar_na(v) for v in joined["node"].tolist()]]
        core_cols = {"parent", "branch.length", "label"}
        extra_cols = ["node"] + [c for c in joined.columns if c not in core_cols and c != "node"]
        td.data = joined._select_cols(extra_cols) if len(extra_cols) > 1 else DataFrame()
        return td

    df_x = x if isinstance(x, DataFrame) else DataFrame(x)
    df_y = y if isinstance(y, DataFrame) else DataFrame(y)
    return _join_dataframes(df_x, df_y, by=by, how="full")


def inner_join(x: Any, y: DataFrame, by: Optional[Union[str, Sequence[str], Dict[str, str]]] = None) -> Any:
    """R dplyr::inner_join."""
    df_x = x if isinstance(x, DataFrame) else DataFrame(x)
    df_y = y if isinstance(y, DataFrame) else DataFrame(y)
    return _join_dataframes(df_x, df_y, by=by, how="inner")


def _join_dataframes(
    df_x: DataFrame,
    df_y: DataFrame,
    by: Optional[Union[str, Sequence[str], Dict[str, str]]] = None,
    how: str = "left",
) -> DataFrame:
    if by is None:
        common = [c for c in df_x.columns if c in df_y.columns]
        if not common:
            if len(df_y.columns) == 0:
                raise ValueError("No columns in right DataFrame for join")
            col0 = df_y.columns[0]
            col0_vals = [v for v in df_y[col0].tolist() if not _is_scalar_na(v)]
            labels_set = (
                {str(v) for v in df_x["label"].tolist() if not _is_scalar_na(v)}
                if "label" in df_x.columns
                else set()
            )
            nodes_set = (
                {v for v in df_x["node"].tolist() if not _is_scalar_na(v)}
                if "node" in df_x.columns
                else set()
            )
            col0_label_hits = sum(1 for v in col0_vals if str(v) in labels_set)
            col0_node_hits = sum(1 for v in col0_vals if v in nodes_set)

            if col0_label_hits > 0 and col0_label_hits >= col0_node_hits:
                by_map = {"label": col0}
            elif col0_node_hits > 0:
                by_map = {"node": col0}
            elif df_y._rownames is not None and len(df_y._rownames) == len(df_y):
                rn_vals = df_y._rownames
                rn_label_hits = sum(1 for v in rn_vals if v in labels_set)
                rn_node_vals = [int(v) if str(v).isdigit() else v for v in rn_vals]
                rn_node_hits = sum(1 for v in rn_node_vals if v in nodes_set)
                if rn_label_hits > 0 and rn_label_hits >= rn_node_hits:
                    df_y = df_y.copy()
                    df_y[".__rowname_key__"] = list(rn_vals)
                    by_map = {"label": ".__rowname_key__"}
                elif rn_node_hits > 0:
                    df_y = df_y.copy()
                    df_y[".__rowname_key__"] = rn_node_vals
                    by_map = {"node": ".__rowname_key__"}
                elif "label" in df_x.columns:
                    by_map = {"label": col0}
                elif "node" in df_x.columns:
                    by_map = {"node": col0}
                else:
                    raise ValueError("No common columns for join")
            elif "label" in df_x.columns:
                by_map = {"label": col0}
            elif "node" in df_x.columns:
                by_map = {"node": col0}
            else:
                raise ValueError("No common columns for join")
        else:
            by_map = {c: c for c in common}
    elif isinstance(by, str):
        by_map = {by: by}
    elif isinstance(by, dict):
        by_map = dict(by)
    else:
        by_map = {c: c for c in by}

    x_keys = list(by_map.keys())
    y_keys = list(by_map.values())

    def make_key(row: Dict[str, Any], keys: List[str]) -> Tuple[Any, ...]:
        parts = []
        for k in keys:
            v = row.get(k)
            if isinstance(v, float) and not math.isnan(v) and v.is_integer():
                v = int(v)
            parts.append(v)
        return tuple(parts)

    x_recs = df_x.to_records()
    y_recs = df_y.to_records()

    y_index: Dict[Tuple[Any, ...], List[int]] = {}
    for j, r in enumerate(y_recs):
        k = make_key(r, y_keys)
        if not any(_is_scalar_na(p) for p in k):
            y_index.setdefault(k, []).append(j)

    out_rows: List[Dict[str, Any]] = []
    matched_y: set = set()
    y_extra_cols = [c for c in df_y.columns if c not in y_keys]

    for rx in x_recs:
        kx = make_key(rx, x_keys)
        j_list = y_index.get(kx, [])
        if j_list:
            for j in j_list:
                matched_y.add(j)
                merged = dict(rx)
                for c in y_extra_cols:
                    val = y_recs[j].get(c)
                    if c not in merged or _is_scalar_na(merged.get(c)) or not _is_scalar_na(val):
                        merged[c] = val
                out_rows.append(merged)
        elif how in ("left", "full"):
            merged = dict(rx)
            for c in y_extra_cols:
                if c not in merged:
                    merged[c] = None
            out_rows.append(merged)

    if how == "full":
        for j, ry in enumerate(y_recs):
            if j not in matched_y:
                merged = {c: None for c in df_x.columns}
                for xk, yk in by_map.items():
                    merged[xk] = ry.get(yk)
                for c in y_extra_cols:
                    merged[c] = ry.get(c)
                out_rows.append(merged)

    cls = type(df_x)
    all_cols = list(df_x.columns) + [c for c in y_extra_cols if c not in df_x.columns]
    res = cls(out_rows, columns=all_cols, attrs=df_x.attrs)
    for k, lvls in df_x._levels.items():
        res._levels[k] = list(lvls)
    for k, lvls in df_y._levels.items():
        if k in res._cols:
            res._levels[k] = list(lvls)
    return res


def bind_rows(*dfs: Any) -> DataFrame:
    """R dplyr::bind_rows equivalent."""
    valid = [ensure_dataframe(d) for d in dfs if d is not None]
    valid = [d for d in valid if len(d) > 0]
    if not valid:
        return DataFrame()
    all_cols: List[str] = []
    for d in valid:
        for c in d.columns:
            if c not in all_cols:
                all_cols.append(c)
    rows: List[Dict[str, Any]] = []
    for d in valid:
        for r in d.to_records():
            rows.append({c: r.get(c, None) for c in all_cols})
    cls = type(valid[0])
    return cls(rows, columns=all_cols, attrs=valid[0].attrs)


def unnest(df: Any, cols: Union[str, Sequence[str]]) -> DataFrame:
    """R tidyr::unnest equivalent."""
    df = ensure_dataframe(df)
    if isinstance(cols, str):
        cols = [cols]
    col_names = [df._resolve_col(c) for c in cols]
    rows: List[Dict[str, Any]] = []
    for r in df.to_records():
        first_val = r.get(col_names[0])
        if isinstance(first_val, (list, tuple, np.ndarray)) and not isinstance(first_val, (str, bytes)):
            m = len(first_val)
            for idx in range(m):
                nr = dict(r)
                for c in col_names:
                    v = r.get(c)
                    nr[c] = v[idx] if isinstance(v, (list, tuple, np.ndarray)) and idx < len(v) else v
                rows.append(nr)
        else:
            rows.append(dict(r))
    cls = type(df)
    return cls(rows, columns=df.columns, attrs=df.attrs)


def read_table(
    file: str,
    header: bool = True,
    sep: Optional[str] = "\t",
    row_names: Optional[Union[int, str]] = None,
    stringsAsFactors: bool = False,
    na_strings: Sequence[str] = ("NA", "", "<NA>"),
    **kwargs: Any,
) -> DataFrame:
    """R's read.table / read.delim / read.csv equivalent."""
    row_names = kwargs.pop("row.names", row_names)
    na_set = set(na_strings)
    with open(file, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\r\n") for line in f if line.strip()]
    if not lines:
        return DataFrame()

    def split_line(ln: str) -> List[str]:
        if sep is None or sep == "":
            return ln.split()
        return [tok.strip('"').strip("'") for tok in ln.split(sep)]

    first_tokens = split_line(lines[0])
    second_tokens = split_line(lines[1]) if len(lines) > 1 else first_tokens

    has_implicit_rownames = header and (len(second_tokens) == len(first_tokens) + 1)
    start_idx = 1 if header else 0

    if header:
        col_names = first_tokens if has_implicit_rownames else first_tokens
    else:
        col_names = [f"V{i+1}" for i in range(len(first_tokens))]

    rnames: List[str] = []
    records: List[Dict[str, Any]] = []

    for ln in lines[start_idx:]:
        toks = split_line(ln)
        if has_implicit_rownames:
            rnames.append(toks[0])
            vals = toks[1:]
        elif row_names is not None and isinstance(row_names, int):
            ridx = row_names - 1 if row_names >= 1 else row_names
            rnames.append(toks[ridx])
            vals = [t for i, t in enumerate(toks) if i != ridx]
            if len(col_names) == len(toks):
                col_names = [c for i, c in enumerate(col_names) if i != ridx]
        else:
            vals = toks

        rec: Dict[str, Any] = {}
        for cname, raw_v in zip(col_names, vals):
            if raw_v in na_set:
                rec[cname] = None
            else:
                try:
                    if "." in raw_v or "e" in raw_v.lower():
                        rec[cname] = float(raw_v)
                    else:
                        rec[cname] = int(raw_v)
                except ValueError:
                    rec[cname] = raw_v
        records.append(rec)

    df = DataFrame(records, columns=col_names)
    if rnames:
        df.rownames = rnames
    return df


def read_csv(file: str, header: bool = True, sep: str = ",", **kwargs: Any) -> DataFrame:
    """R's read.csv equivalent."""
    return read_table(file, header=header, sep=sep, **kwargs)
