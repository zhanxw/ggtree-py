"""
Aesthetic mapping (aes) and R/Python symbolic expression evaluator for ggtree.
Supports:
  - String column names & R expressions: aes(color="group", subset="!isTip", subset='label %in% c("t1", "t3")')
  - Unquoted R-like symbolic variables: aes(color=group, subset=~isTip, subset=(isTip == FALSE), subset=(node == 15))
  - Callables & direct vectors/Series: aes(subset=lambda d: d.node == 15)
"""

from __future__ import annotations

import inspect
import re
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Union
import numpy as np

from .dataframe import ColumnSeries, DataFrame, Factor, FALSE, NA, NULL, TRUE, _is_scalar_na, c, paste, paste0


def _capture_caller_env() -> Dict[str, Any]:
    """Capture caller frame locals/globals (like an R formula/quosure environment)."""
    env: Dict[str, Any] = {}
    frame = inspect.currentframe()
    try:
        while frame is not None:
            fname = frame.f_code.co_filename
            if not (fname.endswith("/aes.py") or fname.endswith("\\aes.py") or fname.endswith("/tidytree.py") or fname.endswith("\\tidytree.py")):
                env.update(frame.f_globals)
                env.update(frame.f_locals)
                break
            frame = frame.f_back
    finally:
        del frame
    return env


class SymExpr:
    """
    Symbolic column/expression proxy allowing unquoted R-like syntax inside aes(...),
    such as:
      aes(color=group)
      aes(label=node, subset=(isTip == False))
      aes(label=node, subset=~isTip)
      aes(subset=(node == 15))
      aes(subset=label.isin(["t1", "t3"]))
    """

    def __init__(self, name: str, fn: Optional[Callable[[DataFrame], Any]] = None):
        self._name = name
        self._fn = fn

    def eval(self, df: DataFrame) -> Any:
        if self._fn is not None:
            return self._fn(df)
        if self._name in df:
            return df[self._name]
        return [self._name] * len(df)

    @property
    def name(self) -> str:
        return self._name

    def isin(self, values: Iterable[Any]) -> "SymExpr":
        vals = list(values)
        return SymExpr(f"{self._name} %in% {vals!r}", lambda df: self.eval(df).isin(vals))

    def _bin(self, other: Any, op_str: str, op_fn: Callable[[Any, Any], Any]) -> "SymExpr":
        def _eval(df: DataFrame) -> Any:
            lhs = self.eval(df)
            rhs = other.eval(df) if isinstance(other, SymExpr) else other
            return op_fn(lhs, rhs)
        other_name = other._name if isinstance(other, SymExpr) else repr(other)
        return SymExpr(f"({self._name} {op_str} {other_name})", _eval)

    def __eq__(self, other: Any) -> "SymExpr":  # type: ignore[override]
        return self._bin(other, "==", lambda a, b: a == b)

    def __ne__(self, other: Any) -> "SymExpr":  # type: ignore[override]
        return self._bin(other, "!=", lambda a, b: a != b)

    def __lt__(self, other: Any) -> "SymExpr":
        return self._bin(other, "<", lambda a, b: a < b)

    def __le__(self, other: Any) -> "SymExpr":
        return self._bin(other, "<=", lambda a, b: a <= b)

    def __gt__(self, other: Any) -> "SymExpr":
        return self._bin(other, ">", lambda a, b: a > b)

    def __ge__(self, other: Any) -> "SymExpr":
        return self._bin(other, ">=", lambda a, b: a >= b)

    def __and__(self, other: Any) -> "SymExpr":
        return self._bin(other, "&", lambda a, b: a & b)

    def __or__(self, other: Any) -> "SymExpr":
        return self._bin(other, "|", lambda a, b: a | b)

    def __invert__(self) -> "SymExpr":
        return SymExpr(f"!{self._name}", lambda df: ~self.eval(df))

    def __add__(self, other: Any) -> "SymExpr":
        return self._bin(other, "+", lambda a, b: a + b)

    def __radd__(self, other: Any) -> "SymExpr":
        return self._bin(other, "+", lambda a, b: b + a)

    def __sub__(self, other: Any) -> "SymExpr":
        return self._bin(other, "-", lambda a, b: a - b)

    def __rsub__(self, other: Any) -> "SymExpr":
        return self._bin(other, "-", lambda a, b: b - a)

    def __mul__(self, other: Any) -> "SymExpr":
        return self._bin(other, "*", lambda a, b: a * b)

    def __rmul__(self, other: Any) -> "SymExpr":
        return self._bin(other, "*", lambda a, b: b * a)

    def __truediv__(self, other: Any) -> "SymExpr":
        return self._bin(other, "/", lambda a, b: a / b)

    def __neg__(self) -> "SymExpr":
        return SymExpr(f"-{self._name}", lambda df: -self.eval(df))

    def __repr__(self) -> str:
        return f"Sym({self._name!r})"


class CallableSymExpr(SymExpr):
    """Symbolic expression that also delegates direct function calls to an underlying function."""

    def __init__(self, name: str, target_fn: Callable[..., Any]):
        super().__init__(name)
        self._target_fn = target_fn

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        return self._target_fn(*args, **kwargs)


class Aes(dict):
    """
    Aesthetic mapping dictionary (compatible with ggplot2::aes).
    Normalizes colour -> color, size/linewidth, etc.
    """

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__()
        self._env: Dict[str, Any] = _capture_caller_env()
        if len(args) >= 1 and "x" not in kwargs:
            self["x"] = args[0]
        if len(args) >= 2 and "y" not in kwargs:
            self["y"] = args[1]
        for k, v in kwargs.items():
            norm_k = "color" if k == "colour" else k
            self[norm_k] = v

    def copy(self) -> "Aes":
        out = Aes()
        out._env = dict(getattr(self, "_env", {}))
        out.update(self)
        return out

    def merge(self, other: Optional[Dict[str, Any]]) -> "Aes":
        out = self.copy()
        if other:
            if isinstance(other, Aes) and getattr(other, "_env", None):
                out._env.update(other._env)
            for k, v in other.items():
                norm_k = "color" if k == "colour" else k
                out[norm_k] = v
        return out


def aes(*args: Any, **kwargs: Any) -> Aes:
    """Construct aesthetic mappings (ggplot2::aes)."""
    out = Aes(*args, **kwargs)
    out._env = _capture_caller_env()
    return out


def modifyList(x: Dict[str, Any], val: Dict[str, Any]) -> Aes:
    """R utils::modifyList equivalent for merging aesthetic mappings."""
    out = x.copy() if isinstance(x, Aes) else Aes(**x)
    if isinstance(val, Aes) and getattr(val, "_env", None):
        out._env.update(val._env)
    for k, v in val.items():
        norm_k = "color" if k == "colour" else k
        out[norm_k] = v
    return out


def _split_top_level(expr: str, sep: str) -> List[str]:
    parts: List[str] = []
    depth = 0
    quote: Optional[str] = None
    last = 0
    i = 0
    n = len(expr)
    k = len(sep)
    while i < n:
        ch = expr[i]
        if quote is not None:
            if ch == "\\" and i + 1 < n:
                i += 2
                continue
            if ch == quote:
                quote = None
        else:
            if ch in ("'", '"'):
                quote = ch
            elif ch in ("(", "[", "{"):
                depth += 1
            elif ch in (")", "]", "}"):
                depth = max(0, depth - 1)
            elif depth == 0 and expr[i : i + k] == sep:
                if sep in ("&", "|") and i + k < n and expr[i + k] == sep:
                    i += k
                    continue
                if sep in ("&", "|") and i > 0 and expr[i - 1] == sep:
                    i += 1
                    continue
                parts.append(expr[last:i].strip())
                i += k
                last = i
                continue
        i += 1
    parts.append(expr[last:].strip())
    return parts


def _parenthesize_logical_ops(expr: str) -> str:
    or_parts = _split_top_level(expr, "|")
    if len(or_parts) > 1:
        return " | ".join(f"({_parenthesize_logical_ops(p)})" for p in or_parts if p)
    and_parts = _split_top_level(expr, "&")
    if len(and_parts) > 1:
        return " & ".join(f"({_parenthesize_logical_ops(p)})" for p in and_parts if p)
    return expr


def _translate_r_expr_to_python(expr: str) -> str:
    """Translate an R filter/aes expression into a valid Python expression over a DataFrame namespace."""
    # Preserve quoted string literals during regex substitutions
    literals: List[str] = []

    def _hide_lit(m: re.Match[str]) -> str:
        literals.append(m.group(0))
        return f"__STR_LIT_{len(literals) - 1}__"

    s = re.sub(r"""(?:"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*')""", _hide_lit, expr.strip())

    # Replace R's && / || with & / |
    s = s.replace("&&", "&").replace("||", "|")
    # Replace R's leading dot identifiers like .id with _dot_id
    s = re.sub(r"(?<![A-Za-z0-9_)\]])\.([A-Za-z_][A-Za-z0-9_]*)", r"_dot_\1", s)
    # Replace R's !var with ~(var) and !( with ~(
    s = re.sub(r"!\s*([A-Za-z_][A-Za-z0-9_.]*)", r"~(\1)", s)
    s = re.sub(r"!\s*\(", r"~(", s)
    # Replace R's `a %in% b` with `(a).isin(b)`
    s = re.sub(
        r"([A-Za-z_][A-Za-z0-9_.]*)\s*%in%\s*(c\([^\)]*\)|\[[^\]]*\]|[A-Za-z_][A-Za-z0-9_.]*)",
        r"(\1).isin(\2)",
        s,
    )
    # Replace R's 1:5 with list(range(1, 5 + 1))
    s = re.sub(r"\b(\d+)\s*:\s*(\d+)\b", r"list(range(\1, \2 + 1))", s)
    # Replace R's is.na(x) with is_na(x)
    s = re.sub(r"\bis\.na\s*\(", "is_na(", s)

    # Restore string literals
    for idx, lit in enumerate(literals):
        s = s.replace(f"__STR_LIT_{idx}__", lit)

    return _parenthesize_logical_ops(s)


def eval_expr(df: DataFrame, expr: Any, extra_env: Optional[Dict[str, Any]] = None) -> Any:
    """
    Evaluate an aesthetic mapping value (SymExpr, string column/expression, callable, or vector)
    against a DataFrame 'df'.
    """
    if expr is None:
        return None
    if isinstance(expr, SymExpr):
        return expr.eval(df)
    if callable(expr):
        return expr(df)
    if isinstance(expr, (ColumnSeries, Factor)):
        return expr
    if isinstance(expr, (list, tuple, np.ndarray)):
        return ColumnSeries(expr)
    if isinstance(expr, str):
        # 1. Direct column name lookup
        if expr in df:
            return df[expr]
        # 2. Evaluate R/Python expression against df columns
        py_expr = _translate_r_expr_to_python(expr)
        env: Dict[str, Any] = {}
        if extra_env:
            env.update(extra_env)
        env.update({
            "c": c,
            "paste": paste,
            "paste0": paste0,
            "TRUE": True,
            "FALSE": False,
            "NA": None,
            "NULL": None,
            "abs": abs,
            "round": np.round,
            "np": np,
        })
        from .dataframe import is_na
        env["is_na"] = is_na
        for col in df.columns:
            env[col] = df[col]
            env[col.replace(".", "_")] = df[col]
            if col.startswith("."):
                env[f"_dot_{col[1:]}"] = df[col]
        try:
            val = eval(py_expr, {"__builtins__": {"range": range, "list": list, "len": len, "abs": abs, "min": min, "max": max, "str": str, "int": int, "float": float, "bool": bool}}, env)
            if isinstance(val, (ColumnSeries, Factor)):
                return val
            if isinstance(val, (list, tuple, np.ndarray)):
                return ColumnSeries(val)
            return ColumnSeries([val] * len(df))
        except Exception:
            # Literal string value
            return ColumnSeries([expr] * len(df))
    return ColumnSeries([expr] * len(df))


def eval_subset(df: DataFrame, subset_expr: Any, extra_env: Optional[Dict[str, Any]] = None) -> Union[ColumnSeries, List[int]]:
    """
    Evaluate `aes(subset=...)` against 'df', supporting:
      - boolean expressions: `isTip == False`, `~isTip`, `"!isTip"`, `"node == 15"`, `'label %in% c("t1", "t3")'`
      - 1-based integer index subsets: `c(1:5)`, `[1, 2, 3, 4, 5]` (as in test-aes-subset.R!)
    """
    res = eval_expr(df, subset_expr, extra_env=extra_env)
    if isinstance(res, ColumnSeries):
        vals = res.tolist()
    elif isinstance(res, (list, tuple, np.ndarray)):
        vals = list(res)
    else:
        vals = [bool(res)] * len(df)

    # Check if 1-based integer row indices (like `aes(subset=c(1:5))` in test-aes-subset.R)
    if len(vals) > 0 and all(isinstance(v, (int, np.integer)) and not isinstance(v, bool) for v in vals):
        # Convert 1-based R row indices to 0-based row indices
        n = len(df)
        return [int(v) - 1 for v in vals if 1 <= int(v) <= n]
    return ColumnSeries([bool(v) if not _is_scalar_na(v) else False for v in vals])


# Pre-defined symbolic variables for unquoted R-style aes(...) usage in Python
isTip = SymExpr("isTip")
node = SymExpr("node")
parent = SymExpr("parent")
label = SymExpr("label")
branch = SymExpr("branch")
branch_length = SymExpr("branch.length")
x = SymExpr("x")
y = SymExpr("y")
angle = SymExpr("angle")
group = SymExpr("group")
rate = SymExpr("rate")
posterior = SymExpr("posterior")
height = SymExpr("height")
value = SymExpr("value")
id = SymExpr("id")
_id = SymExpr(".id")
place = SymExpr("place")
type_ = SymExpr("type")
annote = SymExpr("annote")
image = SymExpr("image")
offset = SymExpr("offset")
ta1 = SymExpr("ta1")
ta2 = SymExpr("ta2")
taxa1 = SymExpr("taxa1")
taxa2 = SymExpr("taxa2")
tree = SymExpr("tree")
v = SymExpr("v")
S = SymExpr("S")
trait = SymExpr("trait")
score = SymExpr("score")
