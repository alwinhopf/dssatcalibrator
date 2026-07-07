"""Write perturbed parameter values into DSSAT input files.

Writers (all column-aware, driven by the ``@`` header so they tolerate DSSAT's
variable-width spacing):

* ``edit_cultivar``  — CROPGRO ``.CUL`` coefficients (the active hemp path).
* ``edit_filex``     — FileX planting details (``PPOP``/``PLRS``/``PLDP``) and
                       initial soil-water multiplier (dormant axes).
* ``edit_soil``      — ``.SOL`` profile: layer-column multipliers + profile
                       scalars (e.g. ``SLPF``); ``extract_soil_profile`` pulls a
                       single profile out of a multi-profile ``SOIL.SOL``.
* ``edit_weather``   — ``.WTH`` daily columns (multiplier / additive offset).
* ``parse_fields``   — read ``ID_SOIL`` / ``WSTA`` from the FileX ``*FIELDS`` row.
"""

from __future__ import annotations

import re
from pathlib import Path


def _fmt(value: float, width: int) -> str:
    """Right-justify a number into ``width`` chars at the highest precision that fits."""
    for dec in range(4, -1, -1):
        s = f"{value:.{dec}f}"
        if len(s) <= width:
            return s.rjust(width)
    return f"{value:.0f}".rjust(width)[:width]


def _fmt_like_token(value: float, token: str) -> str:
    """Format a scalar for a DSSAT cell without needless churn for unchanged cells.

    If the numeric value is unchanged, return the original token byte-for-byte.
    Otherwise, use the highest decimal precision that fits the fixed-width cell.
    For leading-dot source tokens (common in ``.SPE`` files), keep that compact
    style so small decimals retain useful precision in narrow fields.
    """
    width = len(token)
    old = token.strip()

    old_value = _parse(old)
    value = float(value)
    if old_value is not None and abs(value - old_value) <= max(1e-12, 1e-9 * abs(old_value)):
        return token

    leading_dot = bool(re.match(r"^-?\.\d+$", old))
    reserve_leading_blank = bool(token[:1].isspace())
    numeric_width = width - 1 if reserve_leading_blank and width > 1 else width
    for decimals in range(5, -1, -1):
        s = f"{value:.{decimals}f}"
        if leading_dot or (reserve_leading_blank and abs(value) < 1.0):
            if s.startswith("0."):
                s = s[1:]
            elif s.startswith("-0."):
                s = "-." + s.split(".", 1)[1]
        if len(s) <= numeric_width:
            return s.rjust(width)
    return _fmt(value, width)


def _parse(cell: str):
    cell = cell.strip()
    try:
        return float(cell)
    except ValueError:
        return None


def _is_numeric_value(value) -> bool:
    if isinstance(value, bool):
        return False
    try:
        float(value)
    except (TypeError, ValueError):
        return False
    return True


def _fmt_filex_value(value, width: int, old_cell: str = "", *, force_text: bool = False) -> str:
    """Format a FileX replacement cell, supporting numeric values and DSSAT codes."""
    if not force_text and _is_numeric_value(value):
        return _fmt(float(value), width)
    text = str(value).strip()
    if len(text) > width:
        raise ValueError(f"FileX text value '{text}' does not fit in {width} columns.")
    if force_text and old_cell:
        leading = len(old_cell) - len(old_cell.lstrip())
        if leading > 0 and leading + len(text) <= width:
            return (" " * leading) + text + (" " * (width - leading - len(text)))
    align = "left" if old_cell and old_cell.rstrip() == old_cell.strip() else "right"
    return text.ljust(width) if align == "left" else text.rjust(width)


def cultivar_field_map(cul_path: str | Path) -> dict[str, tuple[int, int]]:
    """Return {coefficient_name: (start_col, end_col)} from the ``@VAR#`` header."""
    lines = Path(cul_path).read_text(errors="replace").splitlines()
    header = next((ln for ln in lines if ln.startswith("@VAR#")), None)
    if header is None:
        raise ValueError(f"No '@VAR#' header found in {cul_path}")
    tokens = list(re.finditer(r"\S+", header))
    names = [t.group() for t in tokens]
    eco_i = names.index("ECO#")
    coeff_tokens = tokens[eco_i + 1:]
    boundaries = [tokens[eco_i].end()] + [t.end() for t in coeff_tokens]
    fmap = {}
    for i, t in enumerate(coeff_tokens):
        fmap[t.group()] = (boundaries[i], boundaries[i + 1])
    return fmap


def edit_cultivar(cul_path: str | Path, anchor_code: str, updates: dict[str, float]) -> None:
    """In-place edit of the ``anchor_code`` cultivar row in a CROPGRO ``.CUL``.

    ``updates`` maps coefficient names (e.g. ``CSDL``, ``EM-FL``, ``LFMAX``) to
    new values. Unlisted coefficients are preserved exactly. Only the active
    (non-commented) row whose VAR# equals ``anchor_code`` is modified.
    """
    cul_path = Path(cul_path)
    lines = cul_path.read_text(errors="replace").splitlines()
    fmap = cultivar_field_map(cul_path)
    for name in updates:
        if name not in fmap:
            raise KeyError(f"Coefficient '{name}' not a column in {cul_path.name}; "
                           f"known: {sorted(fmap)}")

    target_idx = next(
        (i for i, ln in enumerate(lines)
         if ln.startswith(anchor_code) and not ln.lstrip().startswith("!")),
        None,
    )
    if target_idx is None:
        raise ValueError(f"Active cultivar row '{anchor_code}' not found in {cul_path.name}")

    line = lines[target_idx]
    last_end = max(hi for _, hi in fmap.values())
    if len(line) < last_end:
        line = line.ljust(last_end)
    chars = list(line)
    for name, (lo, hi) in fmap.items():
        if name in updates:
            old_cell = "".join(chars[lo:hi])
            chars[lo:hi] = list(_fmt_like_token(float(updates[name]), old_cell))
    lines[target_idx] = "".join(chars)
    cul_path.write_text("\n".join(lines) + "\n")


def read_cultivar_values(cul_path: str | Path, anchor_code: str) -> dict[str, float]:
    """Read the current coefficient values for a cultivar row (for verification)."""
    cul_path = Path(cul_path)
    lines = cul_path.read_text(errors="replace").splitlines()
    fmap = cultivar_field_map(cul_path)
    line = next((ln for ln in lines
                 if ln.startswith(anchor_code) and not ln.lstrip().startswith("!")), None)
    if line is None:
        raise ValueError(f"Cultivar '{anchor_code}' not found in {cul_path.name}")
    out = {}
    for name, (lo, hi) in fmap.items():
        out[name] = _parse(line[lo:hi]) if len(line) >= hi else None
    return out


def parse_header_boundaries(header: str) -> dict[str, tuple[int, int]]:
    """Return {token_name: (start_col, end_col)} for columns under a DSSAT header line."""
    tokens = list(re.finditer(r"\S+", header))
    fmap = {}
    for i, t in enumerate(tokens):
        name = t.group().lstrip("@").strip(".")
        # start position is the end of the previous token (or 0)
        start = tokens[i - 1].end() if i > 0 else 0
        end = t.end()
        fmap[name] = (start, end)
    return fmap


def _header_next_token_starts(header: str) -> dict[str, int]:
    """Return the next header-token start column for each normalized field name."""
    tokens = list(re.finditer(r"\S+", header))
    starts = {}
    for i, t in enumerate(tokens):
        name = t.group().lstrip("@").strip(".")
        starts[name] = tokens[i + 1].start() if i + 1 < len(tokens) else t.end()
    return starts


def _is_data_line(line: str) -> bool:
    s = line.strip()
    return bool(s) and not s.startswith(("!", "@", "*")) and bool(re.match(r"[-+]?(?:\d|\.\d)", s))


def _normalize_filex_update(name: str, spec) -> dict:
    if isinstance(spec, dict):
        out = dict(spec)
    else:
        out = {"field": name, "value": spec}
    out.setdefault("field", out.get("dssat", out.get("filex_field", name)))
    out.setdefault("op", "set")
    return out


def _find_filex_section(lines: list[str], section: str) -> tuple[int, int] | None:
    key = section.upper()
    start = next((i for i, ln in enumerate(lines)
                  if ln.lstrip().startswith("*") and key in ln.upper()), None)
    if start is None:
        return None
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("*")), len(lines))
    return start, end


def _apply_filex_section_updates(lines: list[str], updates: list[dict]) -> None:
    for upd in updates:
        section = upd.get("section", "PLANTING DETAILS")
        found = _find_filex_section(lines, section)
        if found is None:
            if upd.get("required", False):
                raise KeyError(f"FileX section '{section}' not found.")
            continue
        start, end = found
        header_prefix = str(upd.get("header_prefix", "") or "")
        field = str(upd["field"])
        raw_value = upd["value"]
        op = str(upd.get("op", "set")).lower()
        force_text = str(upd.get("type", upd.get("format", ""))).lower() in {"text", "str", "string", "raw", "code"}
        row_selector = upd.get("row")
        treatment = upd.get("treatment", upd.get("trt", upd.get("trtno")))

        header_idx = None
        for i in range(start + 1, end):
            stripped = lines[i].lstrip()
            if not stripped.startswith("@"):
                continue
            if header_prefix and not stripped.startswith(header_prefix):
                continue
            fmap = parse_header_boundaries(lines[i])
            if field in fmap:
                header_idx = i
                break
        if header_idx is None:
            if upd.get("required", False):
                raise KeyError(f"FileX field '{field}' not found in section '{section}'.")
            continue

        fmap = parse_header_boundaries(lines[header_idx])
        next_starts = _header_next_token_starts(lines[header_idx])
        lo, hi = fmap[field]
        last_end = max(hi2 for _, hi2 in fmap.values())
        matched = False
        data_row = 0
        for i in range(header_idx + 1, end):
            ln = lines[i]
            if ln.lstrip().startswith("@"):
                break
            if not _is_data_line(ln):
                continue
            data_row += 1
            if row_selector is not None and int(row_selector) != data_row:
                continue
            if treatment is not None and "TRT" in fmap:
                tlo, thi = fmap["TRT"]
                try:
                    if int(float(ln[tlo:thi].strip())) != int(treatment):
                        continue
                except ValueError:
                    continue
            chars = list(ln.ljust(last_end))
            cell_hi = hi
            next_start = next_starts.get(field, hi)
            extended_hi = next_start - 1 if next_start > hi else hi
            if force_text and extended_hi > hi:
                spill = "".join(chars[hi:extended_hi])
                if spill.strip() or len(str(raw_value).strip()) > hi - lo:
                    cell_hi = extended_hi
                    chars = list("".join(chars).ljust(cell_hi))
            old_cell = "".join(chars[lo:cell_hi])
            old = _parse(old_cell)
            if op == "set":
                new = raw_value
            else:
                if old is None:
                    continue
                if not _is_numeric_value(raw_value):
                    raise ValueError(f"FileX operation '{op}' for field '{field}' requires a numeric value.")
                value = float(raw_value)
                if op in {"mult", "multiply"}:
                    new = old * value
                elif op == "add":
                    new = old + value
                else:
                    raise ValueError(f"Unsupported FileX operation '{op}' for field '{field}'.")
            if bool(upd.get("clip_01", False)):
                new = min(max(float(new), 0.0), 1.0)
            chars[lo:cell_hi] = list(_fmt_filex_value(new, cell_hi - lo, old_cell, force_text=force_text))
            lines[i] = "".join(chars)
            matched = True
        if upd.get("required", False) and not matched:
            raise KeyError(f"No FileX data rows matched field '{field}' in section '{section}'.")


def edit_filex(
    filex_path: str | Path,
    mgt_fields: dict[str, float | dict],
    init_updates: dict[str, float | dict],
    section_updates: list[dict] | None = None,
) -> None:
    """Edit management and initial conditions sections in a FileX in-place.

    ``mgt_fields`` keeps the historical shorthand for ``*PLANTING DETAILS``.
    Dict values may target any FileX section, for example
    ``{"irrig_amount": {"section": "IRRIGATION", "field": "IRVAL", "value": 25}}``.
    ``init_updates`` supports the legacy ``initial_soil_water_mult`` plus generic
    ``*INITIAL CONDITIONS`` field updates such as ``SH2O`` or ``SNH4``.
    """
    filex_path = Path(filex_path)
    lines = filex_path.read_text(errors="replace").splitlines()

    generic_updates = list(section_updates or [])
    planting_fields = {}
    for name, spec in (mgt_fields or {}).items():
        upd = _normalize_filex_update(name, spec)
        generic_keys = {
            "section", "header_prefix", "row", "treatment", "trt", "trtno",
            "clip_01", "required", "type", "format",
        }
        use_generic = (
            isinstance(spec, dict)
            and (bool(generic_keys.intersection(upd)) or str(upd.get("op", "set")).lower() != "set")
        )
        if use_generic:
            upd.setdefault("section", "PLANTING DETAILS")
            generic_updates.append(upd)
        else:
            planting_fields[upd["field"]] = upd["value"]

    init_scalar_updates = {}
    for name, spec in (init_updates or {}).items():
        upd = _normalize_filex_update(name, spec)
        if name == "initial_soil_water_mult" and not isinstance(spec, dict):
            init_scalar_updates[name] = spec
            continue
        if isinstance(spec, dict):
            upd.setdefault("section", "INITIAL CONDITIONS")
            if name == "initial_soil_water_mult":
                if not any(k in spec for k in ("field", "dssat", "filex_field")):
                    upd["field"] = "SH2O"
                if "op" not in spec:
                    upd["op"] = "mult"
                upd.setdefault("clip_01", True)
            generic_updates.append(upd)
        else:
            init_scalar_updates[name] = spec

    if generic_updates:
        _apply_filex_section_updates(lines, generic_updates)

    # 1. Edit management (planting details)
    if planting_fields:
        # Find *PLANTING DETAILS section
        planting_sec_idx = None
        for i, ln in enumerate(lines):
            if ln.startswith("*PLANTING DETAILS"):
                planting_sec_idx = i
                break
        
        if planting_sec_idx is not None:
            # Find the @P header line
            header_idx = None
            for i in range(planting_sec_idx + 1, len(lines)):
                if lines[i].startswith("*"):
                    break
                if lines[i].lstrip().startswith("@P"):
                    header_idx = i
                    break
            
            if header_idx is not None:
                fmap = parse_header_boundaries(lines[header_idx])
                # Edit subsequent data rows
                for i in range(header_idx + 1, len(lines)):
                    ln = lines[i]
                    if ln.startswith("*") or ln.lstrip().startswith("@"):
                        break
                    if ln.lstrip().startswith("!") or not ln.strip():
                        continue
                    if ln.strip() and ln.strip()[0].isdigit():
                        # This is a data row, modify the requested fields
                        chars = list(ln)
                        last_end = max(hi for _, hi in fmap.values())
                        if len(chars) < last_end:
                            chars += [" "] * (last_end - len(chars))
                        for name, val in planting_fields.items():
                            if name in fmap:
                                lo, hi = fmap[name]
                                chars[lo:hi] = list(_fmt(float(val), hi - lo))
                        lines[i] = "".join(chars)

    # 2. Edit initial conditions (initial soil water multiplier)
    if init_scalar_updates and "initial_soil_water_mult" in init_scalar_updates:
        mult = float(init_scalar_updates["initial_soil_water_mult"])
        # Find *INITIAL CONDITIONS section
        init_sec_idx = None
        for i, ln in enumerate(lines):
            if ln.startswith("*INITIAL CONDITIONS"):
                init_sec_idx = i
                break
        
        if init_sec_idx is not None:
            # Find the @C header containing SH2O
            header_idx = None
            for i in range(init_sec_idx + 1, len(lines)):
                if lines[i].startswith("*"):
                    break
                if lines[i].lstrip().startswith("@C") and "SH2O" in lines[i]:
                    header_idx = i
                    break
            
            if header_idx is not None:
                fmap = parse_header_boundaries(lines[header_idx])
                if "SH2O" in fmap:
                    lo, hi = fmap["SH2O"]
                    # Edit subsequent data rows
                    for i in range(header_idx + 1, len(lines)):
                        ln = lines[i]
                        if ln.startswith("*") or ln.lstrip().startswith("@"):
                            break
                        if ln.lstrip().startswith("!") or not ln.strip():
                            continue
                        if ln.strip() and ln.strip()[0].isdigit():
                            # Modifying soil water multiplier
                            chars = list(ln)
                            if len(chars) < hi:
                                chars += [" "] * (hi - len(chars))
                            val_str = "".join(chars[lo:hi]).strip()
                            try:
                                val = float(val_str)
                                new_val = val * mult
                                new_val = min(max(new_val, 0.01), 1.0)
                                chars[lo:hi] = list(_fmt(new_val, hi - lo))
                                lines[i] = "".join(chars)
                            except ValueError:
                                pass

    filex_path.write_text("\n".join(lines) + "\n")


def parse_fields(filex_path: str | Path) -> dict:
    """Return field metadata from the FileX ``*FIELDS`` section.

    The FIELDS row mixes right-justified numeric fields with left-justified, often
    over-wide codes (``ID_SOIL`` values are 8 chars under a 7-char header), so we
    align by whitespace-tokenising the header and data rows together rather than
    by fixed column bounds. Coordinates are returned when the optional XCRD/YCRD
    header row is present and not filled with DSSAT's missing-value sentinel.
    """
    def norm_token(token: str) -> str:
        return token.lstrip("@").strip(".")

    def as_float(value):
        try:
            val = float(value)
        except (TypeError, ValueError):
            return None
        return None if val in (-99.0, -999.0) else val

    lines = Path(filex_path).read_text(errors="replace").splitlines()
    try:
        sec = next(i for i, ln in enumerate(lines) if ln.startswith("*FIELDS"))
    except StopIteration:
        return {}

    fields = {}
    i = sec + 1
    while i < len(lines) and not lines[i].startswith("*"):
        if not lines[i].lstrip().startswith("@L"):
            i += 1
            continue
        header = lines[i]
        data = next((lines[j] for j in range(i + 1, len(lines))
                     if not lines[j].lstrip().startswith("@")
                     and not lines[j].startswith("*")
                     and lines[j].strip()
                     and lines[j].strip()[0].isdigit()), None)
        if data is None:
            i += 1
            continue
        htoks = [norm_token(t) for t in header.split()]
        dtoks = data.split()
        if len(dtoks) >= len(htoks):
            fields.update(dict(zip(htoks, dtoks)))
        i += 1

    if not fields:
        return {}

    return {
        "wsta": fields.get("WSTA"),
        "id_soil": fields.get("ID_SOIL"),
        "id_field": fields.get("ID_FIELD"),
        "lat": as_float(fields.get("YCRD")),
        "lon": as_float(fields.get("XCRD")),
        "elev": as_float(fields.get("ELEV")),
    }


def extract_soil_profile(sol_path: str | Path, profile_id: str) -> str:
    """Return the single ``*<profile_id>`` block from a (multi-profile) ``.SOL`` file."""
    lines = Path(sol_path).read_text(errors="replace").splitlines()
    out, capturing = [], False
    for ln in lines:
        if ln.startswith("*"):
            token = ln[1:].split()[0] if len(ln) > 1 and ln[1:].split() else ""
            if capturing:
                break
            if token == profile_id:
                capturing = True
        if capturing:
            out.append(ln)
    while out and (not out[-1].strip() or out[-1].lstrip().startswith("!")):
        out.pop()
    if not out:
        raise ValueError(f"Soil profile '{profile_id}' not found in {Path(sol_path).name}")
    return "\n".join(out) + "\n"


def edit_soil(sol_path: str | Path, profile_id: str,
              layer_mults: dict | None = None, profile_sets: dict | None = None) -> None:
    """Edit one soil profile in a ``.SOL`` in place.

    ``layer_mults`` multiplies layer-table columns across all layers (e.g.
    ``{"SDUL": 1.1, "SLLL": 0.9, "SSAT": 1.05, "SRGF": 1.2}``); ``profile_sets``
    sets profile-level scalars (e.g. ``{"SLPF": 0.9, "SLRO": 70}``). Volumetric
    layer values are clipped to [0, 1]. The profile-level header is the one
    carrying ``SALB``; the layer header is the one carrying ``SDUL``.
    """
    layer_mults = layer_mults or {}
    profile_sets = profile_sets or {}
    sol_path = Path(sol_path)
    lines = sol_path.read_text(errors="replace").splitlines()

    start = next((i for i, ln in enumerate(lines)
                  if ln.startswith("*") and len(ln) > 1 and ln[1:].split()
                  and ln[1:].split()[0] == profile_id), None)
    if start is None:
        raise ValueError(f"Soil profile '{profile_id}' not found in {sol_path.name}")
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("*")), len(lines))

    def _set_row(chars: list, fmap: dict, name: str, text: str):
        lo, hi = fmap[name]
        if len(chars) < hi:
            chars += [" "] * (hi - len(chars))
        chars[lo:hi] = list(text)

    # profile-level scalars (header carrying SALB, not SDUL)
    if profile_sets:
        for j in range(start, end):
            ln = lines[j]
            if ln.lstrip().startswith("@") and "SALB" in ln and "SDUL" not in ln:
                fmap = parse_header_boundaries(ln)
                k = j + 1
                chars = list(lines[k])
                for name, val in profile_sets.items():
                    if name in fmap:
                        lo, hi = fmap[name]
                        _set_row(chars, fmap, name, _fmt(float(val), hi - lo))
                lines[k] = "".join(chars)
                break

    # layer-table multipliers (header carrying SDUL)
    if layer_mults:
        for j in range(start, end):
            ln = lines[j]
            if ln.lstrip().startswith("@") and "SDUL" in ln and "SLB" in ln:
                fmap = parse_header_boundaries(ln)
                for k in range(j + 1, end):
                    row = lines[k]
                    if row.lstrip().startswith(("@", "*")):
                        break
                    if not row.strip() or row.lstrip().startswith("!"):
                        continue
                    if not row.lstrip()[0].isdigit():
                        continue
                    chars = list(row)
                    for name, mult in layer_mults.items():
                        if name not in fmap:
                            continue
                        lo, hi = fmap[name]
                        cur = _parse("".join(chars[lo:hi])) if len(chars) >= hi else None
                        if cur is None or cur == -99:
                            continue
                        new = min(max(cur * float(mult), 0.0), 1.0)
                        _set_row(chars, fmap, name, _fmt(new, hi - lo))
                    lines[k] = "".join(chars)
                break

    sol_path.write_text("\n".join(lines) + "\n")


def edit_weather(wth_path: str | Path, ops: dict) -> None:
    """Edit a ``.WTH`` in place. ``ops`` maps a column to ``(mode, value)`` where
    ``mode`` is ``"mult"`` (value is a factor) or ``"off"`` (value is an additive
    offset), e.g. ``{"SRAD": ("mult", 1.1), "TMAX": ("off", -1.0)}``. Missing
    (-99) cells are left untouched."""
    wth_path = Path(wth_path)
    lines = wth_path.read_text(errors="replace").splitlines()
    hdr_idx = next((i for i, ln in enumerate(lines)
                    if ln.lstrip().startswith("@") and "SRAD" in ln), None)
    if hdr_idx is None:
        return
    fmap = parse_header_boundaries(lines[hdr_idx])

    def _fmt_w(value: float, width: int) -> str:
        """Weather convention: 1 decimal, right-justified; keeps a separator space."""
        s = f"{value:.1f}"
        return s.rjust(width) if len(s) <= width else _fmt(value, width)

    for i in range(hdr_idx + 1, len(lines)):
        ln = lines[i]
        if not ln.strip() or ln.lstrip().startswith(("@", "*", "!", "$")):
            continue
        if not ln.lstrip()[0].isdigit():
            continue
        chars = list(ln)
        for col, (mode, val) in ops.items():
            if col not in fmap:
                continue
            lo, hi = fmap[col]
            cur = _parse("".join(chars[lo:hi])) if len(chars) >= hi else None
            if cur is None or cur == -99:
                continue
            new = cur * float(val) if mode == "mult" else cur + float(val)
            if len(chars) < hi:
                chars += [" "] * (hi - len(chars))
            chars[lo:hi] = list(_fmt_w(new, hi - lo))
        lines[i] = "".join(chars)
    wth_path.write_text("\n".join(lines) + "\n")


def ecotype_field_map(eco_path: str | Path) -> dict[str, tuple[int, int]]:
    """Return {coefficient_name: (start_col, end_col)} from the ``@ECO#`` header."""
    lines = Path(eco_path).read_text(errors="replace").splitlines()
    header = next((ln for ln in lines if ln.startswith("@ECO#")), None)
    if header is None:
        raise ValueError(f"No '@ECO#' header found in {eco_path}")
    tokens = list(re.finditer(r"\S+", header))
    name_i = next(i for i, t in enumerate(tokens) if t.group().startswith("ECONAME"))
    coeff_tokens = tokens[name_i + 1:]
    boundaries = [tokens[name_i].end()] + [t.end() for t in coeff_tokens]
    fmap = {}
    for i, t in enumerate(coeff_tokens):
        fmap[t.group()] = (boundaries[i], boundaries[i + 1])
    return fmap


def edit_ecotype(eco_path: str | Path, anchor_code: str, updates: dict[str, float]) -> None:
    """In-place edit of the ``anchor_code`` ecotype row in a CROPGRO ``.ECO``.

    ``updates`` maps coefficient names to new values. Unlisted coefficients are preserved exactly.
    """
    eco_path = Path(eco_path)
    lines = eco_path.read_text(errors="replace").splitlines()
    fmap = ecotype_field_map(eco_path)
    for name in updates:
        if name not in fmap:
            raise KeyError(f"Coefficient '{name}' not a column in {eco_path.name}; "
                           f"known: {sorted(fmap)}")

    target_idx = next(
        (i for i, ln in enumerate(lines)
         if ln.startswith(anchor_code) and not ln.lstrip().startswith("!")),
        None,
    )
    if target_idx is None:
        raise ValueError(f"Active ecotype row '{anchor_code}' not found in {eco_path.name}")

    line = lines[target_idx]
    last_end = max(hi for _, hi in fmap.values())
    if len(line) < last_end:
        line = line.ljust(last_end)
    chars = list(line)
    for name, (lo, hi) in fmap.items():
        if name in updates:
            old_cell = "".join(chars[lo:hi])
            chars[lo:hi] = list(_fmt_like_token(float(updates[name]), old_cell))
    lines[target_idx] = "".join(chars)
    eco_path.write_text("\n".join(lines) + "\n")


_NUM_RE = re.compile(r"-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def edit_species(spe_path: str | Path, updates: dict[str, float | dict]) -> None:
    """Best-effort in-place edit of named scalar values in a CROPGRO/CERES ``.SPE``.

    ``.SPE`` files are free-form (not column-tabular like ``.CUL``/``.ECO``), so this
    is deliberately conservative: ``updates`` maps a **line key** (a substring that
    uniquely identifies the target line, usually the trailing ``!`` label or the
    leading label token) to a new value, and the **first numeric token** on the
    matched line is replaced in place, preserving width where possible. Species
    coefficients are physiology-defining, so this is gated (``gating.species:
    free``) and intended for new-species adaptation from an analog template only.

    Raises ``KeyError`` if a key matches zero or multiple active lines, so a typo
    never silently edits the wrong constant.
    """
    spe_path = Path(spe_path)
    lines = spe_path.read_text(errors="replace").splitlines()
    for key, val in updates.items():
        token_index = 0
        if isinstance(val, dict):
            token_index = int(val.get("index", val.get("token_index", 0)))
            val = val.get("value")
        matches = [i for i, ln in enumerate(lines)
                   if key in ln and ln.strip() and not ln.lstrip().startswith(("*", "@", "!"))
                   and _NUM_RE.search(ln)]
        if len(matches) != 1:
            raise KeyError(f"Species key '{key}' matched {len(matches)} lines in "
                           f"{spe_path.name} (need exactly 1); refine the key.")
        i = matches[0]
        ln = lines[i]
        nums = list(_NUM_RE.finditer(ln))
        if token_index < 0 or token_index >= len(nums):
            raise KeyError(f"Species key '{key}' token index {token_index} out of range "
                           f"for {spe_path.name} (found {len(nums)} numeric tokens).")
        m = nums[token_index]
        old = m.group(0)
        new = _fmt_like_token(float(val), old) if old else f"{float(val):.3f}"
        lines[i] = ln[:m.start()] + new + ln[m.end():]
    spe_path.write_text("\n".join(lines) + "\n")


def read_cul_calibration_bounds(cul_path: str | Path) -> dict[str, dict[str, float]]:
    """Read DSSAT ``.CUL`` MINIMA/MAXIMA calibration rows into per-coefficient bounds.

    DSSAT cultivar files carry special rows (``VAR#`` ``999991`` = minima, ``999992``
    = maxima) used by GLUE. Returns ``{coeff: {"min": x, "max": y}}`` for the
    coefficients present in both rows. Empty if the file has no such rows.
    """
    cul_path = Path(cul_path)
    lines = cul_path.read_text(errors="replace").splitlines()
    fmap = cultivar_field_map(cul_path)
    rows = {}
    for ln in lines:
        code = ln[:6].strip()
        if code in ("999991", "999992"):
            vals = {}
            for name, (lo, hi) in fmap.items():
                v = _parse(ln[lo:hi]) if len(ln) >= hi else None
                if v is not None:
                    vals[name] = v
            rows["min" if code == "999991" else "max"] = vals
    if "min" not in rows or "max" not in rows:
        return {}
    out = {}
    for name in fmap:
        if name in rows["min"] and name in rows["max"]:
            out[name] = {"min": rows["min"][name], "max": rows["max"][name]}
    return out


def read_ecotype_values(eco_path: str | Path, anchor_code: str) -> dict[str, float]:
    """Read the current coefficient values for an ecotype row (for verification)."""
    eco_path = Path(eco_path)
    lines = eco_path.read_text(errors="replace").splitlines()
    fmap = ecotype_field_map(eco_path)
    line = next((ln for ln in lines
                 if ln.startswith(anchor_code) and not ln.lstrip().startswith("!")), None)
    if line is None:
        raise ValueError(f"Ecotype '{anchor_code}' not found in {eco_path.name}")
    out = {}
    for name, (lo, hi) in fmap.items():
        val = _parse(line[lo:hi])
        if val is not None:
            out[name] = val
    return out


