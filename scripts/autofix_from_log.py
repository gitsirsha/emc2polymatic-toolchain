#!/blue/jsampath/ganguly.sirsha/.conda/envs/polymer-building/bin/python <-- edit this path for python of your environment

################################################################################
#
# autofix_from_log_improper.py
#
# Author: Sirsha Ganguly
#
# Description:
#   Parses a Polymatic/LAMMPS log file for missing bond, angle, dihedral, and
#   improper type errors and automatically patches the LAMMPS data file
#   (data.lmps) by copying force-field parameters from an equivalent existing
#   type. After patching, regenerates types.txt via lmps2types.py so that
#   subsequent Polymatic steps can proceed without manual intervention.
#
#   Designed to run inside a shell loop alongside polym_loop.py:
#
#     while true; do
#       python polym_loop.py 2>&1 | tee polym_loop.log || true
#       if grep -Eq "type '.*' is not defined" polym_loop.log; then
#         python autofix_from_log_improper.py \
#           --pack data.lmps \
#           --log  polym_loop.log \
#           --types_script lmps2types.py
#       else
#         break
#       fi
#     done
#
# Prerequisites / Conventions:
#   1. data.lmps improper-type comments must follow the LAMMPS/topology
#      convention: center atom at index 1  (i.e. "a,CENTER,b,c").
#      EMC writes center-first by default in the parameters section;
#      run swap_improper_comments.py on the initial pack file before 
#      starting the polymerization loop.
#   2. polym.pl must emit "Improper type '...' is not defined" (not
#      "Dihedral type") for improper failures.  See the NOTE in the PR#2
#      description for the required line edits to polym.pl.
#   3. The EQUIV dict below maps linker/reactive atom-type aliases to their
#      base force-field equivalents used only for parameter lookup — the
#      original label is always preserved in the output file.
#
# Changes vs. v1 (day-1 autofix):
#   - Improper errors are now caught directly (no dihedral→improper
#     reclassification needed once polym.pl and data.lmps conventions are
#     aligned).
#   - candidates() for impropers permutes the three non-center atoms while
#     keeping index-1 fixed, matching getImpropType() in Polymatic.pm.
#
################################################################################

import argparse
import re
import shutil
import subprocess
from itertools import permutations

# ---------------------------------------------------------------------------
# Atom-type aliases: used ONLY for parameter lookup, never written to output
# ---------------------------------------------------------------------------
EQUIV = {
    "Lcp": "cp",
    "Lna": "na",
    "c1":  "c",
    "cp0": "cp",
    "c0":  "c",
    "Lc":  "c",
}

PRIMARY = {
    "bond":     "Bond Coeffs",
    "angle":    "Angle Coeffs",
    "dihedral": "Dihedral Coeffs",
    "improper": "Improper Coeffs",
}

GROUPS = {
    "bond":     ["Bond Coeffs"],
    "angle":    ["Angle Coeffs", "BondBond Coeffs", "BondAngle Coeffs"],
    "dihedral": [
        "Dihedral Coeffs",
        "MiddleBondTorsion Coeffs",
        "EndBondTorsion Coeffs",
        "AngleTorsion Coeffs",
        "BondBond13 Coeffs",
        "AngleAngleTorsion Coeffs",
        "BondAngleTorsion Coeffs",
    ],
    "improper": ["Improper Coeffs", "AngleAngle Coeffs"],
}

TYPE_COUNT_KEY = {
    "bond":     "bond",
    "angle":    "angle",
    "dihedral": "dihedral",
    "improper": "improper",
}

ERR_PAT = re.compile(
    r"Error:\s*(Bond|Angle|Dihedral|Improper)\s+type\s+'([^']+)'\s+is not defined",
    re.IGNORECASE,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def canon_tokens(label):
    toks = [t.strip() for t in label.strip().split(",")]
    return [EQUIV.get(t, t) for t in toks]


def canon_label(label):
    return ",".join(canon_tokens(label))


def candidates(label, kind):
    """Return candidate canonical keys to search in the params block."""
    toks = canon_tokens(label)
    keys = [",".join(toks)]
    if kind in ("bond", "angle", "dihedral") and len(toks) in (2, 3, 4):
        rev = ",".join(toks[::-1])
        if rev not in keys:
            keys.append(rev)
    elif kind == "improper" and len(toks) == 4:
        # Center is at index 1 (LAMMPS topology convention).
        # Polymatic.pm::getImpropType keeps index-1 fixed and tries all
        # permutations of the other three atoms — mirror that here.
        center = toks[1]
        rest = [toks[0], toks[2], toks[3]]
        for perm in permutations(rest):
            key = ",".join([perm[0], center, perm[1], perm[2]])
            if key not in keys:
                keys.append(key)
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
    return sum(
        1 for ln in lines[b0:b1]
        if ln.strip() and ln.strip().split()[0].isdigit()
    )


def update_header_type_count(lines, kind, new_count):
    key = TYPE_COUNT_KEY[kind]
    pat = re.compile(r"^\s*(\d+)\s+" + re.escape(key) + r"\s+types\s*$")
    for i, ln in enumerate(lines):
        if pat.match(ln):
            lines[i] = re.sub(r"^\s*\d+", str(new_count), ln)
            return

# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

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
    return {k: v for k, v in out.items() if v}


def patch_pack_inplace(pack_path, missing):
    with open(pack_path, "r") as f:
        lines = f.readlines()

    total_added = 0

    for kind in missing:
        labels = missing[kind]
        primary = PRIMARY[kind]
        sec = find_section(lines, primary)
        if sec is None:
            raise RuntimeError(
                "Section '%s' not found in %s" % (primary, pack_path)
            )
        _, b0, b1 = sec

        canon_to_old = label_map_from_comments(lines, b0, b1)
        tlines_primary = type_lines(lines, b0, b1)
        next_id = (max(tlines_primary.keys()) + 1) if tlines_primary else 1

        existing_comments = set()
        for ln in lines[b0:b1]:
            if "#" in ln:
                existing_comments.add(ln.split("#", 1)[1].strip())

        for miss in labels:
            comment = miss  # write label exactly as reported by polym.pl

            if comment in existing_comments:
                continue

            old_type = None
            for cand in candidates(miss, kind):
                if cand in canon_to_old:
                    old_type = canon_to_old[cand]
                    break
            if old_type is None:
                raise KeyError(
                    "Can't map missing %s '%s' (canonical '%s') to any "
                    "existing parameter entry." % (kind, miss, canon_label(miss))
                )

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
                new_ln = replace_leading_int(tlines[old_type], new_type, comment)
                lines[s1:s1] = [new_ln]

            existing_comments.add(comment)
            total_added += 1

        new_count = count_numeric_lines(lines, primary)
        if new_count:
            update_header_type_count(lines, kind, new_count)

    with open(pack_path, "w") as f:
        f.writelines(lines)

    return total_added


def main():
    ap = argparse.ArgumentParser(
        description="Auto-patch missing LAMMPS type errors from a Polymatic log."
    )
    ap.add_argument("--pack",         default="data.lmps",    help="LAMMPS data file to patch")
    ap.add_argument("--log",          required=True,           help="Polymatic log file to parse")
    ap.add_argument("--types_script", default="lmps2types.py", help="Script to regenerate types.txt")
    args = ap.parse_args()

    with open(args.log, "r") as f:
        log_text = f.read()

    missing = parse_missing(log_text)
    if not missing:
        print("[autofix] No missing-type errors found.")
        return

    bak = args.pack + ".bak"
    try:
        open(bak).close()
    except IOError:
        shutil.copy2(args.pack, bak)
        print("[autofix] Backup created: %s" % bak)

    added = patch_pack_inplace(args.pack, missing)
    print("[autofix] Added %d new type(s) to %s. Kinds patched: %s"
          % (added, args.pack, list(missing.keys())))

    subprocess.check_call(["python", args.types_script, "-i", args.pack])
    print("[autofix] Regenerated types.txt via %s" % args.types_script)


if __name__ == "__main__":
    main()
