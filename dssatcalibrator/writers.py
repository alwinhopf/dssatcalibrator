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


def _parse(cell: str):
    cell = cell.strip()
    try:
        return float(cell)
    except ValueError:
        return None


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
            chars[lo:hi] = list(_fmt(float(updates[name]), hi - lo))
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
        name = t.group()
        if name.startswith("@"):
            name = name[1:]
        # start position is the end of the previous token (or 0)
        start = tokens[i - 1].end() if i > 0 else 0
        end = t.end()
        fmap[name] = (start, end)
    return fmap


def edit_filex(filex_path: str | Path, mgt_fields: dict[str, float], init_updates: dict[str, float]) -> None:
    """Edit management and initial conditions sections in a FileX in-place."""
    filex_path = Path(filex_path)
    lines = filex_path.read_text(errors="replace").splitlines()

    # 1. Edit management (planting details)
    if mgt_fields:
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
                        for name, val in mgt_fields.items():
                            if name in fmap:
                                lo, hi = fmap[name]
                                chars[lo:hi] = list(_fmt(float(val), hi - lo))
                        lines[i] = "".join(chars)

    # 2. Edit initial conditions (initial soil water multiplier)
    if init_updates and "initial_soil_water_mult" in init_updates:
        mult = float(init_updates["initial_soil_water_mult"])
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
    """Return {'wsta', 'id_soil', 'id_field'} from the FileX ``*FIELDS`` header row.

    The FIELDS row mixes right-justified numeric fields with left-justified, often
    over-wide codes (``ID_SOIL`` values are 8 chars under a 7-char header), so we
    align by whitespace-tokenising the header and data rows together rather than
    by fixed column bounds.
    """
    lines = Path(filex_path).read_text(errors="replace").splitlines()
    try:
        sec = next(i for i, ln in enumerate(lines) if ln.startswith("*FIELDS"))
    except StopIteration:
        return {}
    header_idx = next((i for i in range(sec + 1, len(lines))
                       if lines[i].lstrip().startswith("@L") and "ID_SOIL" in lines[i]), None)
    if header_idx is None:
        return {}
    htoks = lines[header_idx].split()
    htoks[0] = htoks[0].lstrip("@")              # "@L" -> "L"
    data = next((lines[i] for i in range(header_idx + 1, len(lines))
                 if lines[i].strip() and lines[i].strip()[0].isdigit()), None)
    if data is None:
        return {}
    dtoks = data.split()
    if len(dtoks) < len(htoks):
        return {}
    m = dict(zip(htoks, dtoks))
    wsta_key = next((k for k in htoks if k.startswith("WSTA")), None)
    return {"wsta": m.get(wsta_key) if wsta_key else None,
            "id_soil": m.get("ID_SOIL"), "id_field": m.get("ID_FIELD")}


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
            chars[lo:hi] = list(_fmt(float(updates[name]), hi - lo))
    lines[target_idx] = "".join(chars)
    eco_path.write_text("\n".join(lines) + "\n")


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



