"""
True vector PDF export and self-contained interactive HTML / Jupyter viewer for ggtree-py.

All functionality works with zero required runtime dependencies beyond the Python
standard library and numpy (using optional cairosvg / rsvg-convert when installed,
and falling back to a pure-Python SVG-to-PDF compiler).
"""

from __future__ import annotations

import html
import json
import math
import re
import shutil
import subprocess
import uuid
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

from .dataframe import _is_scalar_na


# ---------------------------------------------------------------------------
# Color parsing helpers for SVG -> PDF vector compilation
# ---------------------------------------------------------------------------
_PDF_NAMED_COLORS: Dict[str, Tuple[float, float, float]] = {
    "black": (0.0, 0.0, 0.0),
    "white": (1.0, 1.0, 1.0),
    "red": (1.0, 0.0, 0.0),
    "green": (0.0, 0.502, 0.0),
    "blue": (0.0, 0.0, 1.0),
    "yellow": (1.0, 1.0, 0.0),
    "cyan": (0.0, 1.0, 1.0),
    "magenta": (1.0, 0.0, 1.0),
    "grey": (0.502, 0.502, 0.502),
    "gray": (0.502, 0.502, 0.502),
    "darkgrey": (0.663, 0.663, 0.663),
    "darkgray": (0.663, 0.663, 0.663),
    "lightgrey": (0.827, 0.827, 0.827),
    "lightgray": (0.827, 0.827, 0.827),
    "grey50": (0.5, 0.5, 0.5),
    "gray50": (0.5, 0.5, 0.5),
    "grey70": (0.7, 0.7, 0.7),
    "gray70": (0.7, 0.7, 0.7),
    "steelblue": (0.275, 0.510, 0.706),
    "firebrick": (0.698, 0.133, 0.133),
    "darkgreen": (0.0, 0.392, 0.0),
    "darkred": (0.545, 0.0, 0.0),
    "darkblue": (0.0, 0.0, 0.545),
    "orange": (1.0, 0.647, 0.0),
    "purple": (0.502, 0.0, 0.502),
    "pink": (1.0, 0.753, 0.796),
    "brown": (0.647, 0.165, 0.165),
    "navy": (0.0, 0.0, 0.502),
    "lightblue": (0.678, 0.847, 0.902),
    "gold": (1.0, 0.843, 0.0),
    "salmon": (0.980, 0.502, 0.447),
    "tomato": (1.0, 0.388, 0.278),
    "coral": (1.0, 0.498, 0.314),
    "teal": (0.0, 0.502, 0.502),
    "indigo": (0.294, 0.0, 0.510),
    "violet": (0.933, 0.510, 0.933),
    "olive": (0.502, 0.502, 0.0),
    "tan": (0.824, 0.706, 0.549),
    "beige": (0.961, 0.961, 0.863),
    "skyblue": (0.529, 0.808, 0.922),
    "forestgreen": (0.133, 0.545, 0.133),
    "limegreen": (0.196, 0.804, 0.196),
    "darkorange": (1.0, 0.549, 0.0),
    "dodgerblue": (0.118, 0.565, 1.0),
    "royalblue": (0.255, 0.412, 0.882),
    "slategray": (0.439, 0.502, 0.565),
    "slategrey": (0.439, 0.502, 0.565),
}

_GREY_RE = re.compile(r"^gr[ae]y(\d{1,3})$", re.IGNORECASE)
_NUM_RE = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
_URL_ID_RE = re.compile(r"^url\(#([^)]+)\)$")
_PATH_CMD_RE = re.compile(r"([MmLlHhVvCcSsQqTtAaZz])|([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)")


def _parse_svg_color(val: Optional[str]) -> Optional[Tuple[float, float, float]]:
    if not val:
        return None
    s = val.strip().lower()
    if s in ("", "none", "transparent"):
        return None
    if s.startswith("#"):
        hx = s[1:]
        if len(hx) == 3:
            hx = "".join(ch * 2 for ch in hx)
        if len(hx) >= 6:
            try:
                return (
                    int(hx[0:2], 16) / 255.0,
                    int(hx[2:4], 16) / 255.0,
                    int(hx[4:6], 16) / 255.0,
                )
            except ValueError:
                return (0.0, 0.0, 0.0)
    if s.startswith("rgb(") or s.startswith("rgba("):
        nums = [float(m) for m in _NUM_RE.findall(s)]
        if len(nums) >= 3:
            return (
                max(0.0, min(1.0, nums[0] / 255.0)),
                max(0.0, min(1.0, nums[1] / 255.0)),
                max(0.0, min(1.0, nums[2] / 255.0)),
            )
    if s in _PDF_NAMED_COLORS:
        return _PDF_NAMED_COLORS[s]
    m_gr = _GREY_RE.match(s)
    if m_gr:
        v = max(0.0, min(1.0, int(m_gr.group(1)) / 100.0))
        return (v, v, v)
    return (0.13, 0.13, 0.13)


def _parse_float(val: Optional[str], default: float = 0.0) -> float:
    if val is None:
        return default
    s = str(val).strip()
    if not s:
        return default
    s = s.rstrip("px").rstrip("pt").strip()
    m = _NUM_RE.match(s)
    if m:
        try:
            return float(m.group(0))
        except ValueError:
            return default
    return default


def _parse_style_dict(elem: ET.Element) -> Dict[str, str]:
    style_s = elem.attrib.get("style", "")
    out: Dict[str, str] = {}
    for item in style_s.split(";"):
        if ":" in item:
            k, v = item.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def _get_attr(elem: ET.Element, style_map: Dict[str, str], name: str, default: Optional[str] = None) -> Optional[str]:
    if name in elem.attrib:
        return elem.attrib[name]
    if name in style_map:
        return style_map[name]
    return default


def _transform_to_cm(transform_str: str) -> Optional[Tuple[float, float, float, float, float, float]]:
    """Combine SVG transform commands into a single PDF 2D affine matrix (a, b, c, d, e, f)."""
    if not transform_str:
        return None
    cur = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)

    def _mul(m1: Tuple[float, float, float, float, float, float], m2: Tuple[float, float, float, float, float, float]) -> Tuple[float, float, float, float, float, float]:
        a1, b1, c1, d1, e1, f1 = m1
        a2, b2, c2, d2, e2, f2 = m2
        return (
            a1 * a2 + c1 * b2,
            b1 * a2 + d1 * b2,
            a1 * c2 + c1 * d2,
            b1 * c2 + d1 * d2,
            a1 * e2 + c1 * f2 + e1,
            b1 * e2 + d1 * f2 + f1,
        )

    for cmd, args_s in re.findall(r"(translate|scale|rotate|matrix)\s*\(([^)]*)\)", transform_str):
        nums = [float(x) for x in _NUM_RE.findall(args_s)]
        if cmd == "translate" and nums:
            tx = nums[0]
            ty = nums[1] if len(nums) > 1 else 0.0
            cur = _mul(cur, (1.0, 0.0, 0.0, 1.0, tx, ty))
        elif cmd == "scale" and nums:
            sx = nums[0]
            sy = nums[1] if len(nums) > 1 else sx
            cur = _mul(cur, (sx, 0.0, 0.0, sy, 0.0, 0.0))
        elif cmd == "rotate" and nums:
            ang = math.radians(nums[0])
            ca, sa = math.cos(ang), math.sin(ang)
            if len(nums) >= 3:
                cx, cy = nums[1], nums[2]
                cur = _mul(cur, (1.0, 0.0, 0.0, 1.0, cx, cy))
                cur = _mul(cur, (ca, sa, -sa, ca, 0.0, 0.0))
                cur = _mul(cur, (1.0, 0.0, 0.0, 1.0, -cx, -cy))
            else:
                cur = _mul(cur, (ca, sa, -sa, ca, 0.0, 0.0))
        elif cmd == "matrix" and len(nums) >= 6:
            cur = _mul(cur, (nums[0], nums[1], nums[2], nums[3], nums[4], nums[5]))
    return cur


def _svg_arc_to_cubic_segments(
    x0: float,
    y0: float,
    rx: float,
    ry: float,
    phi_deg: float,
    large_arc: float,
    sweep: float,
    x: float,
    y: float,
) -> List[Tuple[float, float, float, float, float, float]]:
    """Convert an SVG elliptical arc (A/a) into PDF cubic Bezier segments."""
    if abs(x0 - x) < 1e-9 and abs(y0 - y) < 1e-9:
        return []
    rx = abs(rx)
    ry = abs(ry)
    if rx < 1e-9 or ry < 1e-9:
        return [(x0, y0, x, y, x, y)]

    phi = math.radians(phi_deg % 360.0)
    cos_p, sin_p = math.cos(phi), math.sin(phi)
    dx2 = (x0 - x) / 2.0
    dy2 = (y0 - y) / 2.0
    x1p = cos_p * dx2 + sin_p * dy2
    y1p = -sin_p * dx2 + cos_p * dy2

    lam = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry)
    if lam > 1.0:
        s_lam = math.sqrt(lam)
        rx *= s_lam
        ry *= s_lam

    rx2, ry2 = rx * rx, ry * ry
    x1p2, y1p2 = x1p * x1p, y1p * y1p
    denom = rx2 * y1p2 + ry2 * x1p2
    if denom <= 1e-12:
        return [(x0, y0, x, y, x, y)]
    num = max(0.0, rx2 * ry2 - denom)
    coef = math.sqrt(num / denom)
    if bool(large_arc) == bool(sweep):
        coef = -coef
    cxp = coef * (rx * y1p) / ry
    cyp = coef * -(ry * x1p) / rx

    cx = cos_p * cxp - sin_p * cyp + (x0 + x) / 2.0
    cy = sin_p * cxp + cos_p * cyp + (y0 + y) / 2.0

    def _vec_angle(ux: float, uy: float, vx: float, vy: float) -> float:
        dot = ux * vx + uy * vy
        nrm = math.hypot(ux, uy) * math.hypot(vx, vy) or 1.0
        ang = math.acos(max(-1.0, min(1.0, dot / nrm)))
        if ux * vy - uy * vx < 0:
            ang = -ang
        return ang

    ux, uy = (x1p - cxp) / rx, (y1p - cyp) / ry
    vx, vy = (-x1p - cxp) / rx, (-y1p - cyp) / ry
    theta1 = _vec_angle(1.0, 0.0, ux, uy)
    dtheta = _vec_angle(ux, uy, vx, vy)
    if not sweep and dtheta > 0:
        dtheta -= 2.0 * math.pi
    elif sweep and dtheta < 0:
        dtheta += 2.0 * math.pi

    n_seg = max(1, int(math.ceil(abs(dtheta) / (math.pi / 2.0))))
    dt = dtheta / n_seg
    segs: List[Tuple[float, float, float, float, float, float]] = []
    for i in range(n_seg):
        t0 = theta1 + i * dt
        t1 = t0 + dt
        alpha = (4.0 / 3.0) * math.tan((t1 - t0) / 4.0)
        c0, s0 = math.cos(t0), math.sin(t0)
        c1, s1 = math.cos(t1), math.sin(t1)
        p1x_e = c0 - alpha * s0
        p1y_e = s0 + alpha * c0
        p2x_e = c1 + alpha * s1
        p2y_e = s1 - alpha * c1
        p3x_e, p3y_e = c1, s1

        def _to_abs(px_e: float, py_e: float) -> Tuple[float, float]:
            sx_r, sy_r = px_e * rx, py_e * ry
            return (
                cx + cos_p * sx_r - sin_p * sy_r,
                cy + sin_p * sx_r + cos_p * sy_r,
            )

        q1 = _to_abs(p1x_e, p1y_e)
        q2 = _to_abs(p2x_e, p2y_e)
        q3 = _to_abs(p3x_e, p3y_e)
        segs.append((q1[0], q1[1], q2[0], q2[1], q3[0], q3[1]))
    return segs


def _emit_svg_path_ops(d_str: str) -> List[str]:
    """Translate an SVG path `d` string into PDF path construction operators."""
    tokens = [m.group(1) or m.group(2) for m in _PATH_CMD_RE.finditer(d_str or "")]
    ops: List[str] = []
    idx = 0
    cmd = ""
    cx, cy = 0.0, 0.0
    sx, sy = 0.0, 0.0
    last_ctrl: Optional[Tuple[float, float]] = None
    last_cmd = ""

    def _next_f() -> float:
        nonlocal idx
        v = float(tokens[idx])
        idx += 1
        return v

    def _is_num(pos: int) -> bool:
        return pos < len(tokens) and not (len(tokens[pos]) == 1 and tokens[pos].isalpha())

    while idx < len(tokens):
        tok = tokens[idx]
        if len(tok) == 1 and tok.isalpha():
            cmd = tok
            idx += 1
            if cmd in ("Z", "z"):
                ops.append("h")
                cx, cy = sx, sy
                last_ctrl = None
                last_cmd = cmd
            continue
        if not cmd:
            break
        if cmd == "M":
            cx, cy = _next_f(), _next_f()
            sx, sy = cx, cy
            ops.append(f"{cx:.3f} {cy:.3f} m")
            cmd = "L"
            last_ctrl = None
        elif cmd == "m":
            cx += _next_f()
            cy += _next_f()
            sx, sy = cx, cy
            ops.append(f"{cx:.3f} {cy:.3f} m")
            cmd = "l"
            last_ctrl = None
        elif cmd == "L":
            cx, cy = _next_f(), _next_f()
            ops.append(f"{cx:.3f} {cy:.3f} l")
            last_ctrl = None
        elif cmd == "l":
            cx += _next_f()
            cy += _next_f()
            ops.append(f"{cx:.3f} {cy:.3f} l")
            last_ctrl = None
        elif cmd == "H":
            cx = _next_f()
            ops.append(f"{cx:.3f} {cy:.3f} l")
            last_ctrl = None
        elif cmd == "h":
            cx += _next_f()
            ops.append(f"{cx:.3f} {cy:.3f} l")
            last_ctrl = None
        elif cmd == "V":
            cy = _next_f()
            ops.append(f"{cx:.3f} {cy:.3f} l")
            last_ctrl = None
        elif cmd == "v":
            cy += _next_f()
            ops.append(f"{cx:.3f} {cy:.3f} l")
            last_ctrl = None
        elif cmd in ("C", "c"):
            rel = cmd == "c"
            x1 = _next_f() + (cx if rel else 0.0)
            y1 = _next_f() + (cy if rel else 0.0)
            x2 = _next_f() + (cx if rel else 0.0)
            y2 = _next_f() + (cy if rel else 0.0)
            x3 = _next_f() + (cx if rel else 0.0)
            y3 = _next_f() + (cy if rel else 0.0)
            ops.append(f"{x1:.3f} {y1:.3f} {x2:.3f} {y2:.3f} {x3:.3f} {y3:.3f} c")
            last_ctrl = (x2, y2)
            cx, cy = x3, y3
        elif cmd in ("S", "s"):
            rel = cmd == "s"
            if last_cmd in ("C", "c", "S", "s") and last_ctrl is not None:
                x1, y1 = 2.0 * cx - last_ctrl[0], 2.0 * cy - last_ctrl[1]
            else:
                x1, y1 = cx, cy
            x2 = _next_f() + (cx if rel else 0.0)
            y2 = _next_f() + (cy if rel else 0.0)
            x3 = _next_f() + (cx if rel else 0.0)
            y3 = _next_f() + (cy if rel else 0.0)
            ops.append(f"{x1:.3f} {y1:.3f} {x2:.3f} {y2:.3f} {x3:.3f} {y3:.3f} c")
            last_ctrl = (x2, y2)
            cx, cy = x3, y3
        elif cmd in ("Q", "q"):
            rel = cmd == "q"
            qx = _next_f() + (cx if rel else 0.0)
            qy = _next_f() + (cy if rel else 0.0)
            x3 = _next_f() + (cx if rel else 0.0)
            y3 = _next_f() + (cy if rel else 0.0)
            c1x = cx + (2.0 / 3.0) * (qx - cx)
            c1y = cy + (2.0 / 3.0) * (qy - cy)
            c2x = x3 + (2.0 / 3.0) * (qx - x3)
            c2y = y3 + (2.0 / 3.0) * (qy - y3)
            ops.append(f"{c1x:.3f} {c1y:.3f} {c2x:.3f} {c2y:.3f} {x3:.3f} {y3:.3f} c")
            last_ctrl = (qx, qy)
            cx, cy = x3, y3
        elif cmd in ("T", "t"):
            rel = cmd == "t"
            if last_cmd in ("Q", "q", "T", "t") and last_ctrl is not None:
                qx, qy = 2.0 * cx - last_ctrl[0], 2.0 * cy - last_ctrl[1]
            else:
                qx, qy = cx, cy
            x3 = _next_f() + (cx if rel else 0.0)
            y3 = _next_f() + (cy if rel else 0.0)
            c1x = cx + (2.0 / 3.0) * (qx - cx)
            c1y = cy + (2.0 / 3.0) * (qy - cy)
            c2x = x3 + (2.0 / 3.0) * (qx - x3)
            c2y = y3 + (2.0 / 3.0) * (qy - y3)
            ops.append(f"{c1x:.3f} {c1y:.3f} {c2x:.3f} {c2y:.3f} {x3:.3f} {y3:.3f} c")
            last_ctrl = (qx, qy)
            cx, cy = x3, y3
        elif cmd in ("A", "a"):
            rel = cmd == "a"
            rx_a = _next_f()
            ry_a = _next_f()
            rot_a = _next_f()
            la_f = _next_f()
            sw_f = _next_f()
            nx = _next_f() + (cx if rel else 0.0)
            ny = _next_f() + (cy if rel else 0.0)
            for seg in _svg_arc_to_cubic_segments(cx, cy, rx_a, ry_a, rot_a, la_f, sw_f, nx, ny):
                ops.append(f"{seg[0]:.3f} {seg[1]:.3f} {seg[2]:.3f} {seg[3]:.3f} {seg[4]:.3f} {seg[5]:.3f} c")
            cx, cy = nx, ny
            last_ctrl = None
        else:
            if _is_num(idx):
                idx += 1
        last_cmd = cmd
    return ops


def _pdf_escape_text(text: str) -> bytes:
    """Encode string as PDF literal bytes in WinAnsi (cp1252) with octal escapes."""
    raw = text.encode("cp1252", errors="replace")
    out = bytearray(b"(")
    for b in raw:
        if b in (40, 41, 92):  # '(', ')', '\'
            out.append(92)
            out.append(b)
        elif b < 32 or b >= 127:
            out.extend(f"\\{b:03o}".encode("ascii"))
        else:
            out.append(b)
    out.append(41)
    return bytes(out)


def svg_to_pdf_pure_python(
    svg: str,
    page_width_pt: Optional[float] = None,
    page_height_pt: Optional[float] = None,
) -> bytes:
    """
    Compile an SVG string directly into a genuine multi-layer vector PDF (%PDF-1.4)
    using only the Python standard library.
    """
    svg_clean = svg.strip()
    if svg_clean.startswith("<?xml"):
        svg_clean = svg_clean.split("?>", 1)[1].strip()
    root = ET.fromstring(svg_clean)

    vb_str = root.attrib.get("viewBox", "")
    vb_nums = [float(x) for x in _NUM_RE.findall(vb_str)]
    raw_w = _parse_float(root.attrib.get("width"), 800.0) or 800.0
    raw_h = _parse_float(root.attrib.get("height"), 600.0) or 600.0
    if len(vb_nums) >= 4 and vb_nums[2] > 0 and vb_nums[3] > 0:
        vb_x, vb_y, vb_w, vb_h = vb_nums[:4]
    else:
        vb_x, vb_y, vb_w, vb_h = 0.0, 0.0, raw_w, raw_h

    pw_pt = float(page_width_pt) if page_width_pt is not None else float(raw_w)
    ph_pt = float(page_height_pt) if page_height_pt is not None else float(raw_h)
    pw_pt = max(pw_pt, 10.0)
    ph_pt = max(ph_pt, 10.0)

    sx = pw_pt / max(vb_w, 1e-6)
    sy = ph_pt / max(vb_h, 1e-6)

    gradients: Dict[str, Dict[str, Any]] = {}
    clip_rects: Dict[str, Tuple[float, float, float, float]] = {}

    for el in root.iter():
        tag = el.tag.split("}", 1)[-1]
        eid = el.attrib.get("id")
        if not eid:
            continue
        if tag in ("linearGradient", "radialGradient"):
            stops: List[Tuple[float, Tuple[float, float, float], float]] = []
            for ch in el:
                if ch.tag.split("}", 1)[-1] != "stop":
                    continue
                st_m = _parse_style_dict(ch)
                off_s = _get_attr(ch, st_m, "offset", "0") or "0"
                off_v = _parse_float(off_s, 0.0)
                if "%" in off_s:
                    off_v /= 100.0
                sc = _parse_svg_color(_get_attr(ch, st_m, "stop-color", "#4682B4")) or (0.27, 0.51, 0.71)
                so = max(0.0, min(1.0, _parse_float(_get_attr(ch, st_m, "stop-opacity", "1"), 1.0)))
                stops.append((max(0.0, min(1.0, off_v)), sc, so))
            stops.sort(key=lambda t: t[0])
            x1g = _parse_float(el.attrib.get("x1"), 0.0)
            y1g = _parse_float(el.attrib.get("y1"), 0.0)
            x2g = _parse_float(el.attrib.get("x2"), 1.0)
            y2g = _parse_float(el.attrib.get("y2"), 0.0)
            gradients[eid] = {"kind": tag, "stops": stops, "dir": (x1g, y1g, x2g, y2g)}
        elif tag == "clipPath":
            for ch in el:
                if ch.tag.split("}", 1)[-1] == "rect":
                    rx_c = _parse_float(ch.attrib.get("x"), 0.0)
                    ry_c = _parse_float(ch.attrib.get("y"), 0.0)
                    rw_c = _parse_float(ch.attrib.get("width"), vb_w)
                    rh_c = _parse_float(ch.attrib.get("height"), vb_h)
                    clip_rects[eid] = (rx_c, ry_c, rw_c, rh_c)
                    break

    extgstate_map: Dict[Tuple[float, float], str] = {}

    def _get_gs_name(ca_fill: float, ca_stroke: float) -> str:
        key = (round(max(0.0, min(1.0, ca_fill)), 2), round(max(0.0, min(1.0, ca_stroke)), 2))
        if key not in extgstate_map:
            extgstate_map[key] = f"GS{len(extgstate_map) + 1}"
        return extgstate_map[key]

    stream_parts: List[bytes] = [
        f"{sx:.6f} 0 0 {-sy:.6f} {-vb_x * sx:.4f} {ph_pt + vb_y * sy:.4f} cm\n".encode("ascii")
    ]

    def _emit_str(s: str) -> None:
        stream_parts.append((s + "\n").encode("ascii"))

    def _apply_clip_and_transform(elem: ET.Element, style_map: Dict[str, str]) -> bool:
        tf_str = elem.attrib.get("transform", "")
        cp_attr = _get_attr(elem, style_map, "clip-path", "") or ""
        cp_m = _URL_ID_RE.match(cp_attr.strip())
        cm = _transform_to_cm(tf_str) if tf_str else None
        cp_box = clip_rects.get(cp_m.group(1)) if cp_m else None
        if cm is None and cp_box is None:
            return False
        _emit_str("q")
        if cm is not None:
            _emit_str(f"{cm[0]:.6f} {cm[1]:.6f} {cm[2]:.6f} {cm[3]:.6f} {cm[4]:.3f} {cm[5]:.3f} cm")
        if cp_box is not None:
            _emit_str(f"{cp_box[0]:.3f} {cp_box[1]:.3f} {cp_box[2]:.3f} {cp_box[3]:.3f} re W n")
        return True

    def _paint_style(
        elem: ET.Element,
        style_map: Dict[str, str],
        inherited: Dict[str, Any],
        default_fill: Optional[str] = "black",
        default_stroke: Optional[str] = "none",
    ) -> Tuple[Optional[Tuple[float, float, float]], Optional[Tuple[float, float, float]], Optional[str], float, float]:
        fill_raw = _get_attr(elem, style_map, "fill", inherited.get("fill", default_fill))
        stroke_raw = _get_attr(elem, style_map, "stroke", inherited.get("stroke", default_stroke))
        fill_op = _parse_float(_get_attr(elem, style_map, "fill-opacity"), 1.0)
        stroke_op = _parse_float(_get_attr(elem, style_map, "stroke-opacity"), 1.0)
        op_all = _parse_float(_get_attr(elem, style_map, "opacity"), 1.0)
        ca_f = fill_op * op_all
        ca_s = stroke_op * op_all

        grad_id = None
        if fill_raw:
            gm = _URL_ID_RE.match(fill_raw.strip())
            if gm:
                grad_id = gm.group(1)
        fill_rgb = None if grad_id else _parse_svg_color(fill_raw)
        stroke_rgb = _parse_svg_color(stroke_raw)
        sw = _parse_float(_get_attr(elem, style_map, "stroke-width", inherited.get("stroke-width", "1")), 1.0)
        if stroke_rgb is None or sw <= 0.0:
            stroke_rgb = None

        if ca_f < 0.999 or ca_s < 0.999:
            gs = _get_gs_name(ca_f, ca_s)
            _emit_str(f"/{gs} gs")
        if fill_rgb is not None:
            _emit_str(f"{fill_rgb[0]:.4f} {fill_rgb[1]:.4f} {fill_rgb[2]:.4f} rg")
        if stroke_rgb is not None:
            _emit_str(f"{stroke_rgb[0]:.4f} {stroke_rgb[1]:.4f} {stroke_rgb[2]:.4f} RG")
            _emit_str(f"{sw:.3f} w")
            lcap = _get_attr(elem, style_map, "stroke-linecap", "")
            cap_code = 1 if lcap == "round" else (2 if lcap == "square" else 0)
            _emit_str(f"{cap_code} J")
            dash_raw = _get_attr(elem, style_map, "stroke-dasharray", "")
            if dash_raw and dash_raw.lower() != "none":
                d_nums = [f"{float(x):.3f}" for x in _NUM_RE.findall(dash_raw)]
                if d_nums:
                    _emit_str(f"[{' '.join(d_nums)}] 0 d")
            else:
                _emit_str("[] 0 d")
        return fill_rgb, stroke_rgb, grad_id, ca_f, ca_s

    def _interp_gradient(stops: List[Tuple[float, Tuple[float, float, float], float]], t: float) -> Tuple[Tuple[float, float, float], float]:
        if not stops:
            return (0.27, 0.51, 0.71), 1.0
        if len(stops) == 1 or t <= stops[0][0]:
            return stops[0][1], stops[0][2]
        if t >= stops[-1][0]:
            return stops[-1][1], stops[-1][2]
        for i in range(len(stops) - 1):
            t0, c0, a0 = stops[i]
            t1, c1, a1 = stops[i + 1]
            if t0 <= t <= t1:
                u = (t - t0) / max(t1 - t0, 1e-9)
                return (
                    (c0[0] + u * (c1[0] - c0[0]), c0[1] + u * (c1[1] - c0[1]), c0[2] + u * (c1[2] - c0[2])),
                    a0 + u * (a1 - a0),
                )
        return stops[-1][1], stops[-1][2]

    def _finish_path(fill_rgb: Optional[Tuple[float, float, float]], stroke_rgb: Optional[Tuple[float, float, float]], close_subpath: bool = False) -> None:
        if fill_rgb is not None and stroke_rgb is not None:
            _emit_str("b*" if close_subpath else "B*")
        elif fill_rgb is not None:
            _emit_str("f*")
        elif stroke_rgb is not None:
            _emit_str("s" if close_subpath else "S")
        else:
            _emit_str("n")

    def _select_pdf_font(weight: str, style_v: str) -> str:
        is_b = "bold" in str(weight).lower()
        is_i = str(style_v).lower() in ("italic", "oblique")
        if is_b and is_i:
            return "F4"
        if is_b:
            return "F2"
        if is_i:
            return "F3"
        return "F1"

    def _walk(elem: ET.Element, inherited: Dict[str, Any]) -> None:
        tag = elem.tag.split("}", 1)[-1]
        cls = elem.attrib.get("class", "")
        if tag in ("defs", "script", "style", "title", "metadata"):
            return
        if any(skip in cls for skip in ("ggtree-tooltip-box", "ggtree-node-overlay")):
            return
        if elem.attrib.get("visibility") == "hidden" or elem.attrib.get("display") == "none":
            return

        style_map = _parse_style_dict(elem)
        pushed = _apply_clip_and_transform(elem, style_map)

        try:
            if tag in ("svg", "g"):
                nxt_inh = dict(inherited)
                for ak in ("fill", "stroke", "stroke-width", "font-family", "font-size", "font-weight", "font-style"):
                    av = _get_attr(elem, style_map, ak)
                    if av is not None:
                        nxt_inh[ak] = av
                for ch in elem:
                    _walk(ch, nxt_inh)
                return

            if tag == "rect":
                rx_p = _parse_float(elem.attrib.get("x"), 0.0)
                ry_p = _parse_float(elem.attrib.get("y"), 0.0)
                w_s = elem.attrib.get("width", "0")
                h_s = elem.attrib.get("height", "0")
                rw_p = vb_w if str(w_s).endswith("%") else _parse_float(w_s, 0.0)
                rh_p = vb_h if str(h_s).endswith("%") else _parse_float(h_s, 0.0)
                if rw_p <= 0 or rh_p <= 0:
                    return
                _emit_str("q")
                fill_rgb, stroke_rgb, grad_id, ca_f, ca_s = _paint_style(elem, style_map, inherited, default_fill="black", default_stroke="none")
                if grad_id and grad_id in gradients and gradients[grad_id]["stops"]:
                    ginfo = gradients[grad_id]
                    stops = ginfo["stops"]
                    x1g, y1g, x2g, y2g = ginfo["dir"]
                    n_steps = 24
                    if abs(y2g - y1g) >= abs(x2g - x1g):
                        # Vertical gradient (e.g. colorbar)
                        rev_v = y1g > y2g
                        step_h = rh_p / n_steps
                        for si in range(n_steps):
                            t_val = (si + 0.5) / n_steps
                            g_t = (1.0 - t_val) if rev_v else t_val
                            gc, ga = _interp_gradient(stops, g_t)
                            if ga * ca_f < 0.999:
                                _emit_str(f"/{_get_gs_name(ga * ca_f, ca_s)} gs")
                            else:
                                _emit_str(f"/{_get_gs_name(1.0, 1.0)} gs")
                            _emit_str(f"{gc[0]:.4f} {gc[1]:.4f} {gc[2]:.4f} rg")
                            _emit_str(f"{rx_p:.3f} {ry_p + si * step_h:.3f} {rw_p:.3f} {step_h + 0.35:.3f} re f")
                    else:
                        rev_h = x1g > x2g
                        step_w = rw_p / n_steps
                        for si in range(n_steps):
                            t_val = (si + 0.5) / n_steps
                            g_t = (1.0 - t_val) if rev_h else t_val
                            gc, ga = _interp_gradient(stops, g_t)
                            if ga * ca_f < 0.999:
                                _emit_str(f"/{_get_gs_name(ga * ca_f, ca_s)} gs")
                            else:
                                _emit_str(f"/{_get_gs_name(1.0, 1.0)} gs")
                            _emit_str(f"{gc[0]:.4f} {gc[1]:.4f} {gc[2]:.4f} rg")
                            _emit_str(f"{rx_p + si * step_w:.3f} {ry_p:.3f} {step_w + 0.35:.3f} {rh_p:.3f} re f")
                    if stroke_rgb is not None:
                        _emit_str(f"/{_get_gs_name(ca_f, ca_s)} gs")
                        _emit_str(f"{rx_p:.3f} {ry_p:.3f} {rw_p:.3f} {rh_p:.3f} re S")
                else:
                    _emit_str(f"{rx_p:.3f} {ry_p:.3f} {rw_p:.3f} {rh_p:.3f} re")
                    _finish_path(fill_rgb, stroke_rgb)
                _emit_str("Q")
                return

            if tag == "line":
                x1_l = _parse_float(elem.attrib.get("x1"), 0.0)
                y1_l = _parse_float(elem.attrib.get("y1"), 0.0)
                x2_l = _parse_float(elem.attrib.get("x2"), 0.0)
                y2_l = _parse_float(elem.attrib.get("y2"), 0.0)
                _emit_str("q")
                _, stroke_rgb, _, _, _ = _paint_style(elem, style_map, inherited, default_fill="none", default_stroke="black")
                if stroke_rgb is not None:
                    _emit_str(f"{x1_l:.3f} {y1_l:.3f} m {x2_l:.3f} {y2_l:.3f} l S")
                _emit_str("Q")
                return

            if tag in ("polyline", "polygon"):
                pts_raw = [float(x) for x in _NUM_RE.findall(elem.attrib.get("points", ""))]
                if len(pts_raw) < 4:
                    return
                _emit_str("q")
                is_poly = tag == "polygon"
                fill_rgb, stroke_rgb, grad_id, _, _ = _paint_style(
                    elem,
                    style_map,
                    inherited,
                    default_fill=("black" if is_poly else "none"),
                    default_stroke="none",
                )
                if grad_id and grad_id in gradients and gradients[grad_id]["stops"]:
                    fill_rgb = gradients[grad_id]["stops"][0][1]
                    _emit_str(f"{fill_rgb[0]:.4f} {fill_rgb[1]:.4f} {fill_rgb[2]:.4f} rg")
                _emit_str(f"{pts_raw[0]:.3f} {pts_raw[1]:.3f} m")
                for pi in range(2, len(pts_raw) - 1, 2):
                    _emit_str(f"{pts_raw[pi]:.3f} {pts_raw[pi + 1]:.3f} l")
                _finish_path(fill_rgb, stroke_rgb, close_subpath=is_poly)
                _emit_str("Q")
                return

            if tag in ("circle", "ellipse"):
                cx_c = _parse_float(elem.attrib.get("cx"), 0.0)
                cy_c = _parse_float(elem.attrib.get("cy"), 0.0)
                if tag == "circle":
                    rx_c = ry_c = _parse_float(elem.attrib.get("r"), 0.0)
                else:
                    rx_c = _parse_float(elem.attrib.get("rx"), 0.0)
                    ry_c = _parse_float(elem.attrib.get("ry"), 0.0)
                if rx_c <= 0 or ry_c <= 0:
                    return
                k = 0.5522847498
                ox, oy = rx_c * k, ry_c * k
                _emit_str("q")
                fill_rgb, stroke_rgb, _, _, _ = _paint_style(elem, style_map, inherited, default_fill="black", default_stroke="none")
                _emit_str(f"{cx_c + rx_c:.3f} {cy_c:.3f} m")
                _emit_str(f"{cx_c + rx_c:.3f} {cy_c + oy:.3f} {cx_c + ox:.3f} {cy_c + ry_c:.3f} {cx_c:.3f} {cy_c + ry_c:.3f} c")
                _emit_str(f"{cx_c - ox:.3f} {cy_c + ry_c:.3f} {cx_c - rx_c:.3f} {cy_c + oy:.3f} {cx_c - rx_c:.3f} {cy_c:.3f} c")
                _emit_str(f"{cx_c - rx_c:.3f} {cy_c - oy:.3f} {cx_c - ox:.3f} {cy_c - ry_c:.3f} {cx_c:.3f} {cy_c - ry_c:.3f} c")
                _emit_str(f"{cx_c + ox:.3f} {cy_c - ry_c:.3f} {cx_c + rx_c:.3f} {cy_c - oy:.3f} {cx_c + rx_c:.3f} {cy_c:.3f} c")
                _finish_path(fill_rgb, stroke_rgb, close_subpath=True)
                _emit_str("Q")
                return

            if tag == "path":
                d_s = elem.attrib.get("d", "")
                p_ops = _emit_svg_path_ops(d_s)
                if not p_ops:
                    return
                _emit_str("q")
                fill_rgb, stroke_rgb, grad_id, _, _ = _paint_style(elem, style_map, inherited, default_fill="black", default_stroke="none")
                if grad_id and grad_id in gradients and gradients[grad_id]["stops"]:
                    fill_rgb = gradients[grad_id]["stops"][0][1]
                    _emit_str(f"{fill_rgb[0]:.4f} {fill_rgb[1]:.4f} {fill_rgb[2]:.4f} rg")
                for op in p_ops:
                    _emit_str(op)
                _finish_path(fill_rgb, stroke_rgb, close_subpath=("h" in p_ops))
                _emit_str("Q")
                return

            if tag == "text":
                tx = _parse_float(elem.attrib.get("x"), 0.0)
                ty = _parse_float(elem.attrib.get("y"), 0.0)
                fsize = _parse_float(_get_attr(elem, style_map, "font-size", inherited.get("font-size", "10")), 10.0)
                anchor = str(_get_attr(elem, style_map, "text-anchor", "start") or "start").lower()
                baseline = str(_get_attr(elem, style_map, "dominant-baseline", "") or "").lower()
                y_adj = (0.35 * fsize) if baseline in ("central", "middle") else 0.0
                if baseline in ("hanging", "text-before-edge"):
                    y_adj = 0.8 * fsize

                # Collect runs preserving tspan baseline-shift / font styling when present
                runs: List[Tuple[str, str, float, float]] = []
                full_text = "".join(elem.itertext())
                if not full_text:
                    return
                has_styled_tspan = any(ch.tag.split("}", 1)[-1] == "tspan" for ch in elem)
                base_weight = str(_get_attr(elem, style_map, "font-weight", inherited.get("font-weight", "normal")) or "normal")
                base_style = str(_get_attr(elem, style_map, "font-style", inherited.get("font-style", "normal")) or "normal")
                if has_styled_tspan:
                    if elem.text:
                        runs.append((elem.text, _select_pdf_font(base_weight, base_style), fsize, 0.0))
                    for ch in elem:
                        if ch.tag.split("}", 1)[-1] == "tspan":
                            ts_m = _parse_style_dict(ch)
                            tw = str(_get_attr(ch, ts_m, "font-weight", base_weight) or base_weight)
                            ts = str(_get_attr(ch, ts_m, "font-style", base_style) or base_style)
                            tfs = _parse_float(_get_attr(ch, ts_m, "font-size"), fsize)
                            bs = str(_get_attr(ch, ts_m, "baseline-shift", "") or "").lower()
                            dy_t = _parse_float(ch.attrib.get("dy"), 0.0)
                            if bs == "super":
                                dy_t -= 0.32 * fsize
                            elif bs == "sub":
                                dy_t += 0.22 * fsize
                            if ch.text:
                                runs.append((ch.text, _select_pdf_font(tw, ts), tfs, dy_t))
                        if ch.tail:
                            runs.append((ch.tail, _select_pdf_font(base_weight, base_style), fsize, 0.0))
                else:
                    runs.append((full_text, _select_pdf_font(base_weight, base_style), fsize, 0.0))

                est_total_w = sum(len(t) * rfs * 0.52 for t, _, rfs, _ in runs)
                dx_anchor = -0.5 * est_total_w if anchor == "middle" else (-est_total_w if anchor == "end" else 0.0)
                _emit_str("q")
                fill_rgb, _, _, _, _ = _paint_style(elem, style_map, inherited, default_fill="black", default_stroke="none")
                if fill_rgb is None:
                    _emit_str("0 0 0 rg")
                cur_x = tx + dx_anchor
                cur_y = ty + y_adj
                _emit_str("BT")
                for t_chunk, f_name, rfs, r_dy in runs:
                    _emit_str(f"/{f_name} {rfs:.2f} Tf")
                    _emit_str(f"1 0 0 -1 {cur_x:.3f} {cur_y + r_dy:.3f} Tm")
                    stream_parts.append(_pdf_escape_text(t_chunk) + b" Tj\n")
                    cur_x += len(t_chunk) * rfs * 0.52
                _emit_str("ET")
                _emit_str("Q")
                return
        finally:
            if pushed:
                _emit_str("Q")

    _walk(root, {"fill": "black", "stroke": "none", "stroke-width": "1"})

    raw_stream = b"".join(stream_parts)
    compressed = zlib.compress(raw_stream, level=6)

    gs_entries = " ".join(
        f"/{name} << /Type /ExtGState /ca {ca_f:.2f} /CA {ca_s:.2f} >>"
        for (ca_f, ca_s), name in extgstate_map.items()
    )
    ext_res = f" /ExtGState << {gs_entries} >>" if gs_entries else ""

    objs: List[bytes] = []
    objs.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    objs.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
    objs.append(
        f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {pw_pt:.2f} {ph_pt:.2f}] "
        f"/Resources 4 0 R /Contents 9 0 R >>".encode("ascii")
    )
    objs.append(
        (
            f"<< /Font << /F1 5 0 R /F2 6 0 R /F3 7 0 R /F4 8 0 R >>"
            f"{ext_res} >>"
        ).encode("ascii")
    )
    for base_font in ("Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique"):
        objs.append(f"<< /Type /Font /Subtype /Type1 /BaseFont /{base_font} /Encoding /WinAnsiEncoding >>".encode("ascii"))
    objs.append(
        b"<< /Length "
        + str(len(compressed)).encode("ascii")
        + b" /Filter /FlateDecode >>\nstream\n"
        + compressed
        + b"\nendstream"
    )

    pdf_out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets: List[int] = [0]
    for obj_id, body in enumerate(objs, start=1):
        offsets.append(len(pdf_out))
        pdf_out.extend(f"{obj_id} 0 obj\n".encode("ascii"))
        pdf_out.extend(body)
        pdf_out.extend(b"\nendobj\n")

    xref_pos = len(pdf_out)
    pdf_out.extend(f"xref\n0 {len(objs) + 1}\n".encode("ascii"))
    pdf_out.extend(b"0000000000 65535 f \n")
    for off in offsets[1:]:
        pdf_out.extend(f"{off:010d} 00000 n \n".encode("ascii"))
    pdf_out.extend(
        f"trailer\n<< /Size {len(objs) + 1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n".encode("ascii")
    )
    return bytes(pdf_out)


def svg_to_pdf_bytes(
    svg: str,
    width: int = 800,
    height: int = 600,
    page_width_pt: Optional[float] = None,
    page_height_pt: Optional[float] = None,
) -> bytes:
    """
    Convert an SVG plot string into true vector PDF bytes.

    Backend priority:
      1. cairosvg (if installed)
      2. rsvg-convert (if installed on PATH)
      3. Pure-Python SVG-to-PDF compiler (always available, zero dependencies)
    """
    pw_pt = float(page_width_pt) if page_width_pt is not None else float(max(int(width), 10))
    ph_pt = float(page_height_pt) if page_height_pt is not None else float(max(int(height), 10))

    try:
        import cairosvg  # type: ignore

        pdf_b = cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), output_width=pw_pt, output_height=ph_pt)
        if isinstance(pdf_b, (bytes, bytearray)) and pdf_b.startswith(b"%PDF-"):
            return bytes(pdf_b)
    except Exception:
        pass

    rsvg_bin = shutil.which("rsvg-convert")
    if rsvg_bin is None:
        for cand in ("/opt/homebrew/bin/rsvg-convert", "/usr/local/bin/rsvg-convert"):
            if Path(cand).exists():
                rsvg_bin = cand
                break
    if rsvg_bin is not None:
        try:
            proc = subprocess.run(
                [
                    rsvg_bin,
                    "-f",
                    "pdf",
                    "-d",
                    "72",
                    "-p",
                    "72",
                    "-w",
                    str(max(int(round(pw_pt)), 1)),
                    "-h",
                    str(max(int(round(ph_pt)), 1)),
                ],
                input=svg.encode("utf-8"),
                capture_output=True,
                check=False,
            )
            if proc.returncode == 0 and proc.stdout.startswith(b"%PDF-"):
                return proc.stdout
        except Exception:
            pass

    return svg_to_pdf_pure_python(svg, page_width_pt=pw_pt, page_height_pt=ph_pt)


# ---------------------------------------------------------------------------
# Interactive SVG node overlay & HTML/Jupyter viewer
# ---------------------------------------------------------------------------
def _clean_meta_value(v: Any) -> Any:
    if v is None or _is_scalar_na(v):
        return None
    if isinstance(v, bool):
        return v
    if isinstance(v, int):
        return v
    if hasattr(v, "item") and callable(v.item):
        try:
            v = v.item()
        except Exception:
            pass
    if isinstance(v, (int, bool)):
        return v
    if isinstance(v, float):
        if not math.isfinite(v):
            return None
        return round(v, 5)
    if isinstance(v, (list, tuple, dict, set)):
        return None
    s = str(v)
    return s if len(s) <= 80 else s[:77] + "..."


def format_tile_bar_tooltip(r: Dict[str, Any], kind: str = "tile") -> str:
    """Build a concise `data-tooltip="..."` attribute for interactive heatmap/bar tiles."""
    parts: List[str] = []
    if kind == "tile":
        lbl = r.get("label") or r.get("id")
        if lbl is not None and not _is_scalar_na(lbl):
            parts.append(str(lbl))
        var = r.get("variable") if not _is_scalar_na(r.get("variable")) else r.get("position")
        if var is not None and not _is_scalar_na(var):
            parts.append(f"col={var}")
        val = r.get("base") if not _is_scalar_na(r.get("base")) else r.get("value", r.get("fill"))
        if val is not None and not _is_scalar_na(val):
            v_str = f"{float(val):.4g}" if isinstance(val, (int, float)) and not isinstance(val, bool) else str(val)
            parts.append(f"value={v_str}")
    else:
        lbl = r.get("label") or r.get("id")
        if lbl is not None and not _is_scalar_na(lbl):
            parts.append(str(lbl))
        pan = r.get("PANEL")
        if pan is not None and not _is_scalar_na(pan) and str(pan) != "Tree":
            parts.append(f"panel={pan}")
        xv = r.get("x")
        if xv is not None and not _is_scalar_na(xv):
            try:
                parts.append(f"x={float(xv):.4g}")
            except Exception:
                pass
        fv = r.get("fill")
        if fv is not None and not _is_scalar_na(fv) and str(fv) != "":
            parts.append(f"fill={fv}")
    if not parts:
        return ""
    return f' data-tooltip="{html.escape(" | ".join(parts), quote=True)}"'


def append_interactive_node_overlay(
    parts: List[str],
    data_df: Any,
    map_pt_fn: Callable[[float, float], Tuple[float, float]],
    layout_str: str = "rectangular",
    is_polar: bool = False,
    polar_arc_d_fn: Optional[Callable[[float, float, float], str]] = None,
) -> None:
    """Append transparent node hit-targets and clade edge overlays for JS interactivity."""
    if data_df is None or len(data_df) == 0:
        return
    cols = getattr(data_df, "columns", ())
    if "node" not in cols or "x" not in cols or "y" not in cols:
        return

    recs = data_df.to_records()
    has_id = ".id" in cols
    node_map: Dict[Any, Dict[str, Any]] = {}
    for r in recs:
        nd = r.get("node")
        if _is_scalar_na(nd) or _is_scalar_na(r.get("x")) or _is_scalar_na(r.get("y")):
            continue
        key = (str(r.get(".id", "")), int(nd)) if has_id else int(nd)
        node_map[key] = r

    if not node_map:
        return

    tree_uid = uuid.uuid4().hex[:8]
    parts.append(f'<g class="ggtree-node-overlay" data-tree-id="{tree_uid}">')

    direct_layouts = {"slanted", "equal_angle", "daylight", "ape", "tree_and_leaf", "unrooted", "external", "radial", "splitnet"}
    for key, r in node_map.items():
        pid_raw = r.get("parent")
        if _is_scalar_na(pid_raw):
            continue
        nid = int(r["node"])
        pid = int(pid_raw)
        if pid == nid:
            continue
        pkey = (str(r.get(".id", "")), pid) if has_id else pid
        prow = node_map.get(pkey)
        if prow is None:
            continue
        rx, ry = float(r["x"]), float(r["y"])
        px0, py0 = float(prow["x"]), float(prow["y"])
        if not all(math.isfinite(v) for v in (rx, ry, px0, py0)):
            continue
        p_par = map_pt_fn(px0, py0)
        p_chi = map_pt_fn(rx, ry)
        if is_polar and polar_arc_d_fn is not None:
            seg_cmds: List[str] = []
            if abs(py0 - ry) > 1e-7:
                try:
                    seg_cmds.append(polar_arc_d_fn(px0, py0, ry))
                except TypeError:
                    seg_cmds.append(polar_arc_d_fn(px0, py0, ry))
            p_cor = map_pt_fn(px0, ry)
            seg_cmds.append(f"M {p_cor[0]:.2f} {p_cor[1]:.2f} L {p_chi[0]:.2f} {p_chi[1]:.2f}")
            d_attr = " ".join(seg_cmds)
        elif layout_str in direct_layouts:
            d_attr = f"M {p_par[0]:.2f} {p_par[1]:.2f} L {p_chi[0]:.2f} {p_chi[1]:.2f}"
        else:
            p_cor = map_pt_fn(px0, ry)
            d_attr = f"M {p_par[0]:.2f} {p_par[1]:.2f} L {p_cor[0]:.2f} {p_cor[1]:.2f} L {p_chi[0]:.2f} {p_chi[1]:.2f}"
        id_suffix = f"{r.get('.id', '')}:" if has_id else ""
        parts.append(
            f'<path class="ggtree-clade-edge" data-node="{id_suffix}{nid}" data-parent="{id_suffix}{pid}" '
            f'd="{d_attr}" fill="none" stroke="none" pointer-events="none"/>'
        )

    hit_r = 4.5 if len(node_map) > 80 else 6.5
    skip_meta = {"x", "y", "branch", "angle", "clade_nodes", "PANEL"}
    for key, r in node_map.items():
        nid = int(r["node"])
        pid = int(r["parent"]) if not _is_scalar_na(r.get("parent")) else nid
        is_tip = bool(r.get("isTip", False))
        lbl_raw = r.get("label")
        lbl_str = "" if (_is_scalar_na(lbl_raw) or lbl_raw is None) else str(lbl_raw)
        rx, ry = float(r["x"]), float(r["y"])
        if not (math.isfinite(rx) and math.isfinite(ry)):
            continue
        px, py = map_pt_fn(rx, ry)
        meta: Dict[str, Any] = {}
        for k, v in r.items():
            if not k or str(k).startswith(("_", ".")) or k in skip_meta:
                continue
            cv = _clean_meta_value(v)
            if cv is not None and cv != "":
                meta[str(k)] = cv
        tip_parts = []
        if lbl_str:
            tip_parts.append(lbl_str)
        tip_parts.append(f"node={nid}")
        if "branch.length" in meta:
            tip_parts.append(f"len={meta['branch.length']}")
        tip_str = " | ".join(tip_parts)
        id_suffix = f"{r.get('.id', '')}:" if has_id else ""
        meta_json = json.dumps(meta, ensure_ascii=False, separators=(",", ":"))
        parts.append(
            f'<circle class="ggtree-hit-node" cx="{px:.2f}" cy="{py:.2f}" r="{hit_r:.1f}" '
            f'fill="transparent" stroke="none" pointer-events="all" style="cursor:pointer" '
            f'data-node="{id_suffix}{nid}" data-parent="{id_suffix}{pid}" data-istip="{1 if is_tip else 0}" '
            f'data-label="{html.escape(lbl_str, quote=True)}" data-tooltip="{html.escape(tip_str, quote=True)}" '
            f'data-meta="{html.escape(meta_json, quote=True)}"/>'
        )

    parts.append("</g>")


def _call_to_svg_interactive(plot: Any, width: int, height: int) -> str:
    try:
        return plot.to_svg(width=width, height=height, _interactive_overlay=True)
    except TypeError:
        return plot.to_svg(width=width, height=height)


def render_interactive_html(
    plot: Any,
    width: int = 860,
    height: int = 620,
    interactive: bool = True,
    full_html: bool = True,
    title: Optional[str] = None,
) -> str:
    """Render a GGTree / PlotList / Composite plot as self-contained interactive HTML."""
    w = max(int(width), 200)
    h = max(int(height), 150)

    if not interactive:
        svg = plot.to_svg(width=w, height=h)
        if not full_html:
            return svg
        return (
            f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
            f'<title>{html.escape(str(title or "ggtree plot"))}</title></head>'
            f'<body style="margin:0;background:#fff;">{svg}</body></html>'
        )

    svg = _call_to_svg_interactive(plot, w, h)
    uid = uuid.uuid4().hex[:8]
    root_id = f"ggtree-iv-{uid}"

    resolved_title = title
    if not resolved_title and hasattr(plot, "labels") and isinstance(plot.labels, dict):
        resolved_title = plot.labels.get("title")
    if not resolved_title:
        resolved_title = "ggtree interactive viewer"
    esc_title = html.escape(str(resolved_title))

    n_tips = 0
    n_nodes = 0
    df = getattr(plot, "data", None)
    if df is not None and hasattr(df, "__len__"):
        n_nodes = len(df)
        if "isTip" in getattr(df, "columns", ()):
            n_tips = sum(1 for v in df["isTip"].tolist() if bool(v))
    badge_txt = f"{n_tips} tips / {n_nodes} nodes" if n_nodes > 0 else "composite figure"

    css = f"""#{root_id}{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;border:1px solid #E5E7EB;border-radius:10px;background:#FFFFFF;color:#111827;overflow:hidden;box-shadow:0 1px 4px rgba(0,0,0,0.06);max-width:100%;user-select:none}}
#{root_id} .iv-bar{{display:flex;flex-wrap:wrap;align-items:center;gap:8px;padding:8px 12px;background:#F8FAFC;border-bottom:1px solid #E5E7EB;font-size:12.5px}}
#{root_id} .iv-title{{font-weight:600;color:#0F172A;margin-right:4px}}
#{root_id} .iv-badge{{background:#E0F2FE;color:#0369A1;border-radius:999px;padding:2px 8px;font-size:11px;font-weight:500}}
#{root_id} .iv-search{{padding:4px 8px;border:1px solid #CBD5E1;border-radius:6px;font-size:12px;width:190px;background:#FFF}}
#{root_id} .iv-search-cnt{{font-size:11.5px;color:#475569;min-width:52px}}
#{root_id} .iv-spacer{{flex:1 1 auto}}
#{root_id} button{{padding:4px 9px;border:1px solid #CBD5E1;border-radius:6px;background:#FFF;color:#1E293B;font-size:12px;cursor:pointer}}
#{root_id} button:hover{{background:#F1F5F9;border-color:#94A3B8}}
#{root_id} .iv-stage{{position:relative;width:100%;height:{h}px;overflow:hidden;background:#FFFFFF;cursor:grab}}
#{root_id} .iv-stage.dragging{{cursor:grabbing}}
#{root_id} .iv-canvas{{width:100%;height:100%;transform-origin:0 0;display:flex;align-items:center;justify-content:center}}
#{root_id} .iv-canvas svg{{max-width:100%;max-height:100%;display:block}}
#{root_id} .iv-tip{{position:fixed;display:none;z-index:99999;pointer-events:none;background:rgba(15,23,42,0.95);color:#F8FAFC;border:1px solid #334155;border-radius:8px;padding:8px 10px;font-size:12px;max-width:300px;box-shadow:0 8px 20px rgba(0,0,0,0.25);line-height:1.35}}
#{root_id} .iv-tip-hdr{{font-weight:600;margin-bottom:4px;display:flex;align-items:center;gap:6px}}
#{root_id} .iv-pill{{font-size:10px;padding:1px 6px;border-radius:999px;background:#2563EB;color:#FFF}}
#{root_id} .iv-pill.int{{background:#7C3AED}}
#{root_id} .iv-kv{{display:grid;grid-template-columns:auto 1fr;gap:2px 10px;font-size:11px}}
#{root_id} .iv-k{{color:#94A3B8}}
#{root_id} .iv-v{{color:#F1F5F9;word-break:break-word}}
#{root_id} .iv-inspector{{padding:7px 12px;background:#F8FAFC;border-top:1px solid #E5E7EB;font-size:12px;color:#334155;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}
#{root_id} .ggtree-hit-node.hovered{{fill:rgba(37,99,235,0.22);stroke:#2563EB;stroke-width:1.6}}
#{root_id} .ggtree-hit-node.clade-selected{{fill:rgba(245,158,11,0.25);stroke:#D97706;stroke-width:1.8}}
#{root_id} .ggtree-hit-node.dimmed{{opacity:0.25}}
#{root_id} .ggtree-hit-node.search-match{{fill:rgba(16,185,129,0.28);stroke:#059669;stroke-width:2.2}}
#{root_id} .ggtree-hit-node.collapsed-root{{fill:rgba(124,58,237,0.25);stroke:#7C3AED;stroke-width:2.2;stroke-dasharray:2,2}}
#{root_id} .ggtree-clade-edge.clade-selected{{stroke:#F59E0B;stroke-width:3;stroke-opacity:0.75}}
#{root_id} .ggtree-clade-edge.collapsed-edge{{stroke:#E2E8F0;stroke-width:3;stroke-opacity:0.85}}"""

    js = f"""(function(){{
  var root = document.getElementById("{root_id}");
  if (!root) return;
  var stage = root.querySelector(".iv-stage");
  var canvas = root.querySelector(".iv-canvas");
  var svg = canvas.querySelector("svg");
  var tipEl = root.querySelector(".iv-tip");
  var inspEl = root.querySelector(".iv-inspector");
  var searchEl = root.querySelector(".iv-search");
  var cntEl = root.querySelector(".iv-search-cnt");
  var scale = 1, tx = 0, ty = 0, dragging = false, moved = false, sx = 0, sy = 0;

  function applyTf(){{
    canvas.style.transform = "translate(" + tx.toFixed(1) + "px," + ty.toFixed(1) + "px) scale(" + scale.toFixed(4) + ")";
  }}
  function zoomBy(factor, cx, cy){{
    var rect = stage.getBoundingClientRect();
    var px = (cx !== undefined ? cx - rect.left : rect.width / 2);
    var py = (cy !== undefined ? cy - rect.top : rect.height / 2);
    var next = Math.max(0.2, Math.min(25, scale * factor));
    var ratio = next / scale;
    tx = px - ratio * (px - tx);
    ty = py - ratio * (py - ty);
    scale = next;
    applyTf();
  }}
  function resetView(){{ scale = 1; tx = 0; ty = 0; applyTf(); }}

  root.querySelector('[data-act="zoom-in"]').onclick = function(){{ zoomBy(1.28); }};
  root.querySelector('[data-act="zoom-out"]').onclick = function(){{ zoomBy(0.78); }};
  root.querySelector('[data-act="zoom-one"]').onclick = function(){{ scale = 1; applyTf(); }};
  root.querySelector('[data-act="reset"]').onclick = function(){{ resetView(); clearSelection(); }};
  root.querySelector('[data-act="clear"]').onclick = function(){{ clearSelection(); searchEl.value = ""; applySearch(); }};

  stage.addEventListener("wheel", function(e){{
    e.preventDefault();
    zoomBy(e.deltaY < 0 ? 1.14 : 0.88, e.clientX, e.clientY);
  }}, {{passive: false}});

  stage.addEventListener("mousedown", function(e){{
    if (e.button !== 0) return;
    dragging = true; moved = false;
    sx = e.clientX - tx; sy = e.clientY - ty;
    stage.classList.add("dragging");
  }});
  window.addEventListener("mousemove", function(e){{
    if (!dragging) return;
    var nx = e.clientX - sx, ny = e.clientY - sy;
    if (Math.hypot(nx - tx, ny - ty) > 3) moved = true;
    tx = nx; ty = ny; applyTf();
  }});
  window.addEventListener("mouseup", function(){{
    dragging = false;
    stage.classList.remove("dragging");
  }});

  function treeDescendants(overlay, startId, includeSelf){{
    var childrenMap = {{}};
    overlay.querySelectorAll(".ggtree-hit-node").forEach(function(n){{
      var nid = n.getAttribute("data-node"), pid = n.getAttribute("data-parent");
      if (!childrenMap[pid]) childrenMap[pid] = [];
      if (nid !== pid) childrenMap[pid].push(nid);
    }});
    var visited = new Set();
    var q = [startId];
    if (includeSelf) visited.add(startId);
    while (q.length){{
      var cur = q.shift();
      var ch = childrenMap[cur] || [];
      for (var i = 0; i < ch.length; i++){{
        if (!visited.has(ch[i])){{
          visited.add(ch[i]);
          q.push(ch[i]);
        }}
      }}
    }}
    return visited;
  }}

  function clearSelection(){{
    root.querySelectorAll(".ggtree-hit-node").forEach(function(n){{
      n.classList.remove("clade-selected", "dimmed");
    }});
    root.querySelectorAll(".ggtree-clade-edge").forEach(function(ed){{
      ed.classList.remove("clade-selected");
    }});
    inspEl.textContent = "Hover nodes/tiles for details. Click a node to highlight its clade, double-click internal nodes to toggle collapse, scroll to zoom, drag to pan.";
  }}

  function selectClade(hitNode){{
    var overlay = hitNode.closest(".ggtree-node-overlay");
    if (!overlay) return;
    var nid = hitNode.getAttribute("data-node");
    var desc = treeDescendants(overlay, nid, true);
    var tipCnt = 0;
    overlay.querySelectorAll(".ggtree-hit-node").forEach(function(n){{
      if (desc.has(n.getAttribute("data-node"))){{
        n.classList.add("clade-selected");
        n.classList.remove("dimmed");
        if (n.getAttribute("data-istip") === "1") tipCnt++;
      }} else {{
        n.classList.remove("clade-selected");
        n.classList.add("dimmed");
      }}
    }});
    overlay.querySelectorAll(".ggtree-clade-edge").forEach(function(ed){{
      ed.classList.toggle("clade-selected", desc.has(ed.getAttribute("data-node")));
    }});
    var lbl = hitNode.getAttribute("data-label") || ("Node " + nid);
    var metaTxt = "";
    try {{
      var m = JSON.parse(hitNode.getAttribute("data-meta") || "{{}}");
      metaTxt = Object.keys(m).slice(0, 6).map(function(k){{ return k + "=" + m[k]; }}).join(" | ");
    }} catch (err) {{}}
    inspEl.textContent = "Selected: " + lbl + " (" + desc.size + " nodes, " + tipCnt + " tips)" + (metaTxt ? " — " + metaTxt : "");
  }}

  var collapsedRoots = new Set();
  function toggleCollapse(hitNode){{
    if (hitNode.getAttribute("data-istip") === "1") return;
    var overlay = hitNode.closest(".ggtree-node-overlay");
    if (!overlay) return;
    var tid = overlay.getAttribute("data-tree-id") || "0";
    var nid = hitNode.getAttribute("data-node");
    var key = tid + "::" + nid;
    if (collapsedRoots.has(key)) collapsedRoots.delete(key);
    else collapsedRoots.add(key);

    var hiddenSet = new Set();
    var rootKeys = new Set();
    collapsedRoots.forEach(function(ck){{
      var parts = ck.split("::");
      if (parts[0] === tid){{
        rootKeys.add(parts[1]);
        treeDescendants(overlay, parts[1], false).forEach(function(d){{ hiddenSet.add(d); }});
      }}
    }});
    overlay.querySelectorAll(".ggtree-hit-node").forEach(function(n){{
      var curId = n.getAttribute("data-node");
      n.classList.toggle("collapsed-root", rootKeys.has(curId));
      n.style.display = hiddenSet.has(curId) ? "none" : "";
    }});
    overlay.querySelectorAll(".ggtree-clade-edge").forEach(function(ed){{
      ed.classList.toggle("collapsed-edge", hiddenSet.has(ed.getAttribute("data-node")));
    }});
  }}

  function showTooltip(e, target){{
    if (dragging) return;
    var isNode = target.classList && target.classList.contains("ggtree-hit-node");
    if (isNode){{
      target.classList.add("hovered");
      var lbl = target.getAttribute("data-label") || ("Node " + target.getAttribute("data-node"));
      var isTip = target.getAttribute("data-istip") === "1";
      var htmlOut = '<div class="iv-tip-hdr"><span>' + escapeHtml(lbl) + '</span><span class="iv-pill ' + (isTip ? "" : "int") + '">' + (isTip ? "Tip" : "Internal") + '</span></div><div class="iv-kv">';
      try {{
        var meta = JSON.parse(target.getAttribute("data-meta") || "{{}}");
        var keys = Object.keys(meta);
        for (var i = 0; i < Math.min(keys.length, 10); i++){{
          htmlOut += '<span class="iv-k">' + escapeHtml(keys[i]) + '</span><span class="iv-v">' + escapeHtml(String(meta[keys[i]])) + '</span>';
        }}
      }} catch (err) {{}}
      htmlOut += '</div>';
      tipEl.innerHTML = htmlOut;
    }} else {{
      var msg = target.getAttribute("data-tooltip") || "";
      if (!msg) return;
      tipEl.innerHTML = '<div class="iv-v">' + escapeHtml(msg) + '</div>';
    }}
    tipEl.style.display = "block";
    var txp = Math.min(window.innerWidth - 315, e.clientX + 14);
    var typ = Math.max(8, Math.min(window.innerHeight - 160, e.clientY + 14));
    tipEl.style.left = txp + "px";
    tipEl.style.top = typ + "px";
  }}
  function escapeHtml(s){{
    return String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;");
  }}

  var lastHovered = null;
  stage.addEventListener("mousemove", function(e){{
    if (lastHovered){{ lastHovered.classList.remove("hovered"); lastHovered = null; }}
    var t = e.target && e.target.closest ? e.target.closest(".ggtree-hit-node, [data-tooltip]") : null;
    if (!t){{ tipEl.style.display = "none"; return; }}
    if (t.classList && t.classList.contains("ggtree-hit-node")) lastHovered = t;
    showTooltip(e, t);
  }});
  stage.addEventListener("mouseleave", function(){{
    if (lastHovered){{ lastHovered.classList.remove("hovered"); lastHovered = null; }}
    tipEl.style.display = "none";
  }});
  stage.addEventListener("click", function(e){{
    if (moved) return;
    var hn = e.target && e.target.closest ? e.target.closest(".ggtree-hit-node") : null;
    if (hn) selectClade(hn);
    else clearSelection();
  }});
  stage.addEventListener("dblclick", function(e){{
    var hn = e.target && e.target.closest ? e.target.closest(".ggtree-hit-node") : null;
    if (hn){{ e.preventDefault(); toggleCollapse(hn); }}
  }});

  function applySearch(){{
    var q = (searchEl.value || "").trim().toLowerCase();
    var nodes = root.querySelectorAll(".ggtree-hit-node");
    if (!q){{
      nodes.forEach(function(n){{ n.classList.remove("search-match"); }});
      cntEl.textContent = "";
      return;
    }}
    var hits = 0;
    nodes.forEach(function(n){{
      var hay = ((n.getAttribute("data-label") || "") + " " + (n.getAttribute("data-meta") || "")).toLowerCase();
      var ok = hay.indexOf(q) !== -1;
      n.classList.toggle("search-match", ok);
      if (ok) hits++;
    }});
    cntEl.textContent = hits + (hits === 1 ? " match" : " matches");
  }}
  searchEl.addEventListener("input", applySearch);

  root.querySelector('[data-act="dl-svg"]').onclick = function(){{
    if (!svg) return;
    var clone = svg.cloneNode(true);
    var ov = clone.querySelectorAll(".ggtree-node-overlay, .ggtree-tooltip-box");
    ov.forEach(function(el){{ el.parentNode && el.parentNode.removeChild(el); }});
    var ser = new XMLSerializer().serializeToString(clone);
    var blob = new Blob([ser], {{type: "image/svg+xml;charset=utf-8"}});
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "ggtree_plot.svg";
    a.click();
  }};

  root.querySelector('[data-act="dl-png"]').onclick = function(){{
    if (!svg) return;
    var clone = svg.cloneNode(true);
    clone.querySelectorAll(".ggtree-node-overlay, .ggtree-tooltip-box").forEach(function(el){{ el.parentNode && el.parentNode.removeChild(el); }});
    var ser = new XMLSerializer().serializeToString(clone);
    var url = URL.createObjectURL(new Blob([ser], {{type: "image/svg+xml;charset=utf-8"}}));
    var img = new Image();
    img.onload = function(){{
      var c = document.createElement("canvas");
      c.width = (img.width || {w}) * 2;
      c.height = (img.height || {h}) * 2;
      var ctx = c.getContext("2d");
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, c.width, c.height);
      ctx.drawImage(img, 0, 0, c.width, c.height);
      URL.revokeObjectURL(url);
      c.toBlob(function(b){{
        var a = document.createElement("a");
        a.href = URL.createObjectURL(b);
        a.download = "ggtree_plot.png";
        a.click();
      }}, "image/png");
    }};
    img.src = url;
  }};
}})();"""

    widget_html = (
        f'<style>{css}</style>'
        f'<div id="{root_id}" class="ggtree-interactive-widget">'
        f'<div class="iv-bar">'
        f'<span class="iv-title">{esc_title}</span>'
        f'<span class="iv-badge">{html.escape(badge_txt)}</span>'
        f'<input class="iv-search" type="search" placeholder="Search taxa, node ID, or metadata...">'
        f'<span class="iv-search-cnt"></span>'
        f'<span class="iv-spacer"></span>'
        f'<button type="button" data-act="zoom-out">−</button>'
        f'<button type="button" data-act="zoom-one">100%</button>'
        f'<button type="button" data-act="zoom-in">+</button>'
        f'<button type="button" data-act="reset">Reset</button>'
        f'<button type="button" data-act="clear">Clear Selection</button>'
        f'<button type="button" data-act="dl-svg">Download SVG</button>'
        f'<button type="button" data-act="dl-png">Download PNG</button>'
        f'</div>'
        f'<div class="iv-stage"><div class="iv-canvas">{svg}</div></div>'
        f'<div class="iv-inspector">Hover nodes/tiles for details. Click a node to highlight its clade, double-click internal nodes to toggle collapse, scroll to zoom, drag to pan.</div>'
        f'<div class="iv-tip"></div>'
        f'</div>'
        f'<script>{js}</script>'
    )

    if not full_html:
        return widget_html
    return (
        f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
        f'<meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>{esc_title}</title></head>'
        f'<body style="margin:12px;background:#F1F5F9;">{widget_html}</body></html>'
    )


class InteractivePlot:
    """
    Jupyter-ready interactive plot container with pan/zoom, hover tooltips,
    search, clade highlighting, and vector SVG/PDF + PNG export.
    """

    def __init__(
        self,
        plot: Any,
        width: int = 860,
        height: int = 620,
        title: Optional[str] = None,
    ):
        self.plot = plot
        self.width = int(width)
        self.height = int(height)
        self.title = title

    def to_html(self, full_html: bool = True) -> str:
        return render_interactive_html(
            self.plot,
            width=self.width,
            height=self.height,
            interactive=True,
            full_html=full_html,
            title=self.title,
        )

    def _repr_html_(self) -> str:
        doc = self.to_html(full_html=True)
        return (
            f'<iframe srcdoc="{html.escape(doc, quote=True)}" width="100%" '
            f'height="{self.height + 96}" '
            f'style="border:1px solid #E5E7EB;border-radius:10px;background:#fff;display:block;" '
            f'sandbox="allow-scripts allow-downloads"></iframe>'
        )

    def to_svg(self, width: Optional[int] = None, height: Optional[int] = None) -> str:
        w = int(width) if width is not None else self.width
        h = int(height) if height is not None else self.height
        return self.plot.to_svg(width=w, height=h)

    def to_png(self, width: Optional[int] = None, height: Optional[int] = None) -> bytes:
        w = int(width) if width is not None else self.width
        h = int(height) if height is not None else self.height
        return self.plot.to_png(width=w, height=h)

    def to_pdf(
        self,
        width: Optional[int] = None,
        height: Optional[int] = None,
        page_width_pt: Optional[float] = None,
        page_height_pt: Optional[float] = None,
    ) -> bytes:
        w = int(width) if width is not None else self.width
        h = int(height) if height is not None else self.height
        if hasattr(self.plot, "to_pdf"):
            try:
                return self.plot.to_pdf(width=w, height=h, page_width_pt=page_width_pt, page_height_pt=page_height_pt)
            except TypeError:
                return self.plot.to_pdf(width=w, height=h)
        return svg_to_pdf_bytes(self.to_svg(width=w, height=h), width=w, height=h, page_width_pt=page_width_pt, page_height_pt=page_height_pt)

    def save(self, filename: Union[str, Path], width: float = 8.6, height: float = 6.2, dpi: int = 100, **kwargs: Any) -> str:
        from .plot import ggsave

        return ggsave(filename, plot=self.plot, width=width, height=height, dpi=dpi, interactive=True, **kwargs)

    def __repr__(self) -> str:
        return f"<InteractivePlot width={self.width} height={self.height} plot={self.plot!r}>"


def girafe(
    ggobj: Any = None,
    width_svg: float = 8.6,
    height_svg: float = 6.2,
    dpi: int = 100,
    title: Optional[str] = None,
    **kwargs: Any,
) -> InteractivePlot:
    """
    Create an interactive Jupyter/HTML figure from a ggtree plot (ggiraph::girafe-compatible).
    """
    if ggobj is None:
        from .plot import last_plot

        ggobj = last_plot()
    if ggobj is None:
        raise ValueError("No plot provided and no last_plot() available.")
    w_px = int(round(width_svg * dpi)) if float(width_svg) <= 40.0 else int(round(width_svg))
    h_px = int(round(height_svg * dpi)) if float(height_svg) <= 40.0 else int(round(height_svg))
    return InteractivePlot(ggobj, width=w_px, height=h_px, title=title)
