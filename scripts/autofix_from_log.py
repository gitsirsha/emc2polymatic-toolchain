#!/usr/bin/env python
# Python 2.7 compatible

import argparse
import re
import shutil
import subprocess

EQUIV = {"Lc1": "c1","Lcp": "cp","cp0": "cp","cp1": "cp","cp2": "cp","cp3": "cp",}

PRIMARY = {
    "bond": "Bond Coeffs",
    "angle": "Angle Coeffs",
    "dihedral": "Dihedral Coeffs",
    "improper": "Improper Coeffs",
}

GROUPS = {
    "bond": ["Bond Coeffs"],
    "angle": ["Angle Coeffs", "BondBond Coeffs", "BondAngle Coeffs"],
    "dihedral": [
        "Dihedral Coeffs",
        "MiddleBondTorsion Coeffs",
        "EndBondTorsion Coeffs",
        "AngleTorsion Coeffs",
        "BondBond13 Coeffs",
        "AngleAngleTorsion Coeffs",
        "BondAngleTorsion Coeffs",
    ],
    "improper": ["Improper Coeffs"],
}

TYPE_COUNT_KEY = {"bond": "bond", "angle": "angle", "dihedral": "dihedral", "improper": "improper"}

ERR_PAT = re.compile(
    r"Error:\s*(Bond|Angle|Dihedral|Improper)\s+type\s+'([^']+)'\s+is not defined",
    re.IGNORECASE
)

def canon_tokens(label):
    toks = [t.strip() for t in label.strip().split(",")]
    out = []
    for t in toks:
        out.append(EQUIV.get(t, t))
    return out

def canon_label(label):
    return ",".join(canon_tokens(label))

def candidates(label, kind):
    toks = canon_tokens(label)
    keys = [",".join(toks)]
    if kind == "bond" and len(toks) == 2:
        keys.append(",".join(toks[::-1]))
    elif kind == "angle" and len(toks) == 3:
        keys.append(",".join(toks[::-1]))
    elif kind == "dihedral" and len(toks) == 4:
        keys.append(",".join(toks[::-1]))
    return keys

def find_section(lines, header):
    for i, ln in enumerate(lines):
        if ln.strip() == header:
            j = i + 1
            while j < len(lines) and lines[j].strip() == "":
                j += 1
            k = j
            while k < len(lines) and lines[k].strip() != "":
                k += 1
            return (i, j, k)
    return None

def type_lines(lines, b0, b1):
    out = {}
    for ln in lines[b0:b1]:
        s = ln.strip()
        if not s:
            continue
        first = s.split()[0]
        if first.isdigit():
            out[int(first)] = ln
    return out

def label_map_from_comments(lines, b0, b1):
    out = {}
    for ln in lines[b0:b1]:
        if "#" not in ln:
            continue
        left, right = ln.split("#", 1)
        right = right.strip()
        left = left.strip()
        if not right or not left:
            continue
        first = left.split()[0]
        if not first.isdigit():
            continue
        out[canon_label(right)] = int(first)
    return out

def replace_leading_int(line, new_id, new_comment):
    orig = line.rstrip("\n")
    m = re.match(r"^(\s*)\d+(\s+.*)$", orig)
    if not m:
        raise ValueError("Bad coeff line (no leading int): %r" % line)
    indent = m.group(1)
    rest = m.group(2)

    if "#" in rest:
        before_hash, _ = rest.split("#", 1)
        return "%s%d%s# %s\n" % (indent, new_id, before_hash, new_comment)
    else:
        return "%s%d%s # %s\n" % (indent, new_id, rest, new_comment)

def count_numeric_lines(lines, header):
    sec = find_section(lines, header)
    if sec is None:
        return 0
    _, b0, b1 = sec
    n = 0
    for ln in lines[b0:b1]:
        parts = ln.strip().split()
        if parts and parts[0].isdigit():
            n += 1
    return n

def update_header_type_count(lines, kind, new_count):
    key = TYPE_COUNT_KEY[kind]
    pat = re.compile(r"^\s*(\d+)\s+" + re.escape(key) + r"\s+types\s*$")
    for i, ln in enumerate(lines):
        if pat.match(ln):
            lines[i] = re.sub(r"^\s*\d+", str(new_count), ln)
            return

def parse_missing(log_text):
    out = {"bond": [], "angle": [], "dihedral": [], "improper": []}
    seen = set()
    for m in ERR_PAT.finditer(log_text):
        kind = m.group(1).lower()
        label = m.group(2).strip()
        key = (kind, label)
        if key not in seen:
            seen.add(key)
            out[kind].append(label)
    # drop empties
    out2 = {}
    for k in out:
        if out[k]:
            out2[k] = out[k]
    return out2

def patch_pack_inplace(pack_path, missing):
    with open(pack_path, "r") as f:
        lines = f.readlines()

    total_added = 0

    for kind in missing:
        labels = missing[kind]
        primary = PRIMARY[kind]
        sec = find_section(lines, primary)
        if sec is None:
            raise RuntimeError("Section '%s' not found in %s" % (primary, pack_path))
        _, b0, b1 = sec

        canon_to_old = label_map_from_comments(lines, b0, b1)
        tlines_primary = type_lines(lines, b0, b1)
        next_id = (max(tlines_primary.keys()) + 1) if tlines_primary else 1

        existing_comments = set()
        for ln in lines[b0:b1]:
            if "#" in ln:
                existing_comments.add(ln.split("#", 1)[1].strip())

        for miss in labels:
            if miss in existing_comments:
                continue

            old_type = None
            for cand in candidates(miss, kind):
                if cand in canon_to_old:
                    old_type = canon_to_old[cand]
                    break
            if old_type is None:
                raise KeyError("Can't map missing %s '%s' (canonical '%s')" %
                               (kind, miss, canon_label(miss)))

            new_type = next_id
            next_id += 1

            for header in GROUPS[kind]:
                sec2 = find_section(lines, header)
                if sec2 is None:
                    continue
                _, s0, s1 = sec2
                tlines = type_lines(lines, s0, s1)
                if old_type not in tlines:
                    continue
                new_ln = replace_leading_int(tlines[old_type], new_type, miss)
                lines[s1:s1] = [new_ln]

            existing_comments.add(miss)
            total_added += 1

        new_count = count_numeric_lines(lines, primary)
        if new_count:
            update_header_type_count(lines, kind, new_count)

    with open(pack_path, "w") as f:
        f.writelines(lines)

    return total_added

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", default="data.lmps")
    ap.add_argument("--log", required=True)
    ap.add_argument("--types_script", default="lmps2types.py")
    args = ap.parse_args()

    pack = args.pack
    log_path = args.log

    with open(log_path, "r") as f:
        log_text = f.read()

    missing = parse_missing(log_text)
    if not missing:
        print("[autofix] No missing-type errors found.")
        return

    bak = pack + ".bak"
    try:
        with open(bak, "r"):
            pass
    except IOError:
        shutil.copy2(pack, bak)

    added = patch_pack_inplace(pack, missing)
    print("[autofix] Added %d new types into %s. Missing was: %s" % (added, pack, missing))

    # regenerate types.txt using your script
    subprocess.check_call(["python", args.types_script, "-i", pack])
    print("[autofix] Regenerated types.txt")

if __name__ == "__main__":
    main()
