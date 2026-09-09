#!/usr/bin/env python

################################################################################
#
# fix_permuted_coeffs.py
#
# Author: Sirsha Ganguly
#
# Description:
#   Companion to autofix_from_log.py. Reads autofix_permuted.log
#   and fixes cross-term coefficients in data.lmps for any type that was
#   added via a permutation or reversal match (i.e. "not as is").
#
#   When autofix copies a parameter entry via a permutation/reversal, the
#   main coeffs section is always safe (symmetric), but cross-term sections
#   are directional and must be corrected:
#
#   Angle reversal (a,b,c → c,b,a):
#     BondBond Coeffs    [M r1 r2]          → r1 <-> r2
#     BondAngle Coeffs   [N1 N2 r1 r2]      → N1<->N2, r1<->r2
#
#   Dihedral reversal (a,b,c,d → d,c,b,a):
#     EndBondTorsion Coeffs  [B1 B2 B3 C1 C2 C3 r1 r3] → B↔C, r1↔r3
#     AngleTorsion Coeffs    [D1 D2 D3 E1 E2 E3 θ1 θ2] → D↔E, θ1↔θ2
#     AngleAngleTorsion Coeffs [M θ1 θ2]               → θ1↔θ2
#     BondBond13 Coeffs      [N r1 r3]                  → r1↔r3
#
#   Improper permutation (center J fixed, I/K/L permuted):
#     AngleAngle Coeffs  [M1 M2 M3 θ1 θ2 θ3]
#       → full remap derived from the specific permutation applied
#         (computed analytically, not from a lookup table)
#
#   Safe sections (no swap needed):
#     Angle Coeffs, Dihedral Coeffs, Improper Coeffs,
#     MiddleBondTorsion Coeffs, BondBond Coeffs (M is symmetric)
#     — MBT uses middle bond r_jk which is invariant under reversal,
#       and cos(nφ) is even so A1/A2/A3 are unaffected.
#
# Usage:
#   python fix_permuted_coeffs.py \
#       --pack   data.lmps \
#       --log    autofix_permuted.log
#
# Run this once after autofix_from_log.py + polym_loop has
# finished (or after every autofix pass if you prefer). It is idempotent
# provided the log is not re-appended between runs — use --clear_log to
# truncate autofix_permuted.log after a successful fix so entries are not
# applied twice.
#
################################################################################

import argparse
import ast
import os
import re
import shutil
from itertools import permutations

# ---------------------------------------------------------------------------
# Atom-type aliases — must match EQUIV in autofix_from_log.py
# ---------------------------------------------------------------------------
EQUIV = {
    "Lcp": "cp",
    "Lna": "na",
    "c1":  "c",
    "cp0": "cp",
    "c0":  "c",
    "Lc":  "c",
}  # change this according to your system

# ---------------------------------------------------------------------------
# Section helpers (shared logic with autofix_from_log.py)
# ---------------------------------------------------------------------------

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


def get_coeff_line_idx(lines, header, type_id):
    """Return the line index in `lines` for the given type_id within header."""
    sec = find_section(lines, header)
    if sec is None:
        return None
    _, b0, b1 = sec
    for i in range(b0, b1):
        s = lines[i].strip()
        if not s:
            continue
        parts = s.split()
        if parts[0].isdigit() and int(parts[0]) == type_id:
            return i
    return None


def swap_cols(line, i, j):
    """
    Swap value columns i and j (1-indexed, after the leading type int).
    Preserves the comment if present.
    """
    orig = line.rstrip("\n")
    comment = ""
    if "#" in orig:
        orig, comment = orig.split("#", 1)
        comment = " #" + comment

    parts = orig.split()
    # parts[0] is type id, parts[1:] are value columns (1-indexed)
    # i and j are 1-indexed into the value columns
    parts[i], parts[j] = parts[j], parts[i]
    return " ".join(parts) + comment + "\n"


def swap_col_blocks(line, slice_a, slice_b):
    """
    Swap two equal-length blocks of value columns (1-indexed slices).
    e.g. slice_a=(1,4), slice_b=(4,7) swaps cols 1-3 with cols 4-6.
    """
    orig = line.rstrip("\n")
    comment = ""
    if "#" in orig:
        orig, comment = orig.split("#", 1)
        comment = " #" + comment

    parts = orig.split()
    # parts[0] is type id; value cols start at index 1
    a_start, a_end = slice_a   # 1-indexed, exclusive end
    b_start, b_end = slice_b
    block_a = parts[a_start:a_end]
    block_b = parts[b_start:b_end]
    parts[a_start:a_end] = block_b
    parts[b_start:b_end] = block_a
    return " ".join(parts) + comment + "\n"


def remap_cols(line, remap):
    """
    Reorder value columns (1-indexed) according to remap list.
    remap[i] = j means new column i comes from old column j (both 1-indexed).
    """
    orig = line.rstrip("\n")
    comment = ""
    if "#" in orig:
        orig, comment = orig.split("#", 1)
        comment = " #" + comment

    parts = orig.split()
    vals = parts[1:]  # strip type id
    new_vals = [vals[j - 1] for j in remap]
    return parts[0] + " " + " ".join(new_vals) + comment + "\n"


# ---------------------------------------------------------------------------
# Improper AngleAngle remap
# ---------------------------------------------------------------------------

def improper_aa_remap(miss_label, matched_label):
    """
    Compute the column remap for AngleAngle Coeffs when an improper type
    was matched via a permutation of the non-center atoms.

    Convention (LAMMPS class2, center = atom J = index 1 in quadruplet):
      quadruplet order: I, J, K, L
      θ1 = angle(I-J-K)
      θ2 = angle(I-J-L)
      θ3 = angle(K-J-L)
      M1 couples θ1 × θ3
      M2 couples θ1 × θ2
      M3 couples θ2 × θ3

    The data file line (after type id) has columns:
      [M1, M2, M3, θ1, θ2, θ3]  →  1-indexed: 1 2 3 4 5 6

    Bug fixes vs naive implementation:
      1. EQUIV is applied to miss_label tokens before comparison so that
         Polymatic aliases (e.g. Lna→na, c1→c) don't cause a KeyError
         when matched against the canonical form stored in matched_label.
      2. Position indices (0,1,2 for I,K,L) are used instead of atom type
         name frozensets so that degenerate cases (two non-center atoms of
         the same type, e.g. hc,c,hc,cp) never cause a dict key collision.

    Returns a 6-element remap list (1-indexed) for use with remap_cols().
    """
    # Bug fix 1: apply EQUIV to miss_label so both sides are canonical.
    new_canon = [EQUIV.get(t.strip(), t.strip())
                 for t in miss_label.split(",")]
    old_canon = [EQUIV.get(t.strip(), t.strip())
                 for t in matched_label.split(",")]

    # Non-center positions (J = index 1 is center, fixed):
    # new_nc[i] and old_nc[i] are the canonical atom types at positions
    # 0=I, 1=K, 2=L in the respective quadruplets.
    new_nc = [new_canon[0], new_canon[2], new_canon[3]]
    old_nc = [old_canon[0], old_canon[2], old_canon[3]]

    # Bug fix 2: find the permutation of POSITION INDICES [0,1,2] that maps
    # new non-center positions onto old non-center positions.
    # pos_map[j] = i means old position j came from new position i.
    # Using indices avoids any frozenset collision when two non-center atoms
    # share the same type (e.g. two hc atoms).
    pos_map = None
    for perm in permutations([0, 1, 2]):
        if all(new_nc[perm[j]] == old_nc[j] for j in range(3)):
            pos_map = list(perm)
            break
    if pos_map is None:
        raise ValueError(
            "Cannot find position permutation mapping '%s' → '%s'. "
            "Check EQUIV dict." % (miss_label, matched_label)
        )

    # Old θ angles in terms of old position-index pairs:
    #   θ1_old = angle(I_old-J-K_old) = positions (0,1)
    #   θ2_old = angle(I_old-J-L_old) = positions (0,2)
    #   θ3_old = angle(K_old-J-L_old) = positions (1,2)
    # Express as frozensets of the CORRESPONDING NEW position indices
    # (via pos_map) so lookups are in a consistent space.
    old_theta = {
        frozenset([pos_map[0], pos_map[1]]): 0,   # θ1_old
        frozenset([pos_map[0], pos_map[2]]): 1,   # θ2_old
        frozenset([pos_map[1], pos_map[2]]): 2,   # θ3_old
    }

    # New θ angles in terms of new position-index pairs:
    #   θ1_new = (0,1), θ2_new = (0,2), θ3_new = (1,2)
    new_theta_pairs = [
        frozenset([0, 1]),
        frozenset([0, 2]),
        frozenset([1, 2]),
    ]

    # theta_map[i] = j means new θ_{i+1} = old θ_{j+1}  (0-based)
    theta_map = [old_theta[pair] for pair in new_theta_pairs]

    # M couplings are identified by frozensets of θ indices (0-based).
    old_m = {
        frozenset([0, 2]): 0,   # M1_old couples θ1×θ3
        frozenset([0, 1]): 1,   # M2_old couples θ1×θ2
        frozenset([1, 2]): 2,   # M3_old couples θ2×θ3
    }

    new_m_pairs = [
        frozenset([theta_map[0], theta_map[2]]),   # new M1 couples θ1_new×θ3_new
        frozenset([theta_map[0], theta_map[1]]),   # new M2 couples θ1_new×θ2_new
        frozenset([theta_map[1], theta_map[2]]),   # new M3 couples θ2_new×θ3_new
    ]

    m_map = [old_m[pair] for pair in new_m_pairs]

    # Build 1-indexed remap for columns [M1,M2,M3,θ1,θ2,θ3]:
    # columns 1..3 = M1..M3, columns 4..6 = θ1..θ3
    remap = (
        [m + 1 for m in m_map    ] +
        [t + 4 for t in theta_map]
    )
    return remap


# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------

LOG_PAT = re.compile(
    r"kind=(\w+)\s+new_type=(\d+)\s+cand_idx=(\d+)\s+label=('.*?')\s+matched=('.*?')\s*$"
)


def parse_permuted_log(log_path):
    """
    Returns list of dicts:
      {kind, new_type, cand_idx, label, matched}
    """
    entries = []
    seen = set()
    if not os.path.exists(log_path):
        return entries
    with open(log_path, "r") as f:
        for line in f:
            m = LOG_PAT.search(line)
            if not m:
                continue
            kind      = m.group(1)
            new_type  = int(m.group(2))
            cand_idx  = int(m.group(3))
            label     = ast.literal_eval(m.group(4))
            matched   = ast.literal_eval(m.group(5))
            key = (kind, new_type)
            if key not in seen:
                seen.add(key)
                entries.append(dict(
                    kind=kind, new_type=new_type, cand_idx=cand_idx,
                    label=label, matched=matched
                ))
    return entries


# ---------------------------------------------------------------------------
# Per-kind fixers
# ---------------------------------------------------------------------------

def fix_angle(lines, new_type, **_):
    """Reversal a,b,c → c,b,a: swap directional cross-term columns."""
    fixes = {
        # BondBond Coeffs: type | M | r1 | r2
        #   r1 = i-j bond, r2 = j-k bond → swap on reversal
        "BondBond Coeffs":  [(2, 3)],
        # BondAngle Coeffs: type | N1 | N2 | r1 | r2
        #   N1 couples r_ij side, N2 couples r_jk side → swap N1↔N2 and r1↔r2
        "BondAngle Coeffs": [(1, 2), (3, 4)],
    }
    for header, swaps in fixes.items():
        idx = get_coeff_line_idx(lines, header, new_type)
        if idx is None:
            continue
        for (ci, cj) in swaps:
            lines[idx] = swap_cols(lines[idx], ci, cj)


def fix_dihedral(lines, new_type, **_):
    """Reversal a,b,c,d → d,c,b,a: swap directional cross-term columns."""
    # EndBondTorsion Coeffs: type | B1 B2 B3 | C1 C2 C3 | r1 | r3
    #   B's couple to r_ij end, C's couple to r_kl end → swap blocks + r1↔r3
    idx = get_coeff_line_idx(lines, "EndBondTorsion Coeffs", new_type)
    if idx is not None:
        lines[idx] = swap_col_blocks(lines[idx], (1, 4), (4, 7))
        lines[idx] = swap_cols(lines[idx], 7, 8)

    # AngleTorsion Coeffs: type | D1 D2 D3 | E1 E2 E3 | θ1 | θ2
    #   D's couple to θ_ijk, E's to θ_jkl → swap blocks + θ1↔θ2
    idx = get_coeff_line_idx(lines, "AngleTorsion Coeffs", new_type)
    if idx is not None:
        lines[idx] = swap_col_blocks(lines[idx], (1, 4), (4, 7))
        lines[idx] = swap_cols(lines[idx], 7, 8)

    # AngleAngleTorsion Coeffs: type | M | θ1 | θ2
    #   M is fine, but θ1 and θ2 reference the two end angles → swap
    idx = get_coeff_line_idx(lines, "AngleAngleTorsion Coeffs", new_type)
    if idx is not None:
        lines[idx] = swap_cols(lines[idx], 2, 3)

    # BondBond13 Coeffs: type | N | r1 | r3
    #   r1 = i-j bond, r3 = k-l bond → swap on reversal
    idx = get_coeff_line_idx(lines, "BondBond13 Coeffs", new_type)
    if idx is not None:
        lines[idx] = swap_cols(lines[idx], 2, 3)


def fix_improper(lines, new_type, label, matched, **_):
    """
    Permutation of non-center atoms: remap AngleAngle Coeffs columns.
    The remap is computed analytically from the atom ordering in label
    vs matched — no lookup table, handles all 6 permutations correctly.
    """
    idx = get_coeff_line_idx(lines, "AngleAngle Coeffs", new_type)
    if idx is None:
        return
    remap = improper_aa_remap(label, matched)
    lines[idx] = remap_cols(lines[idx], remap)


FIXERS = {
    "angle":    fix_angle,
    "dihedral": fix_dihedral,
    "improper": fix_improper,
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Fix cross-term coefficients in data.lmps for types that were "
            "added via permutation/reversal match by autofix_from_log.py."
        )
    )
    ap.add_argument("--pack",      default="data.lmps",           help="LAMMPS data file to patch")
    ap.add_argument("--log",       default="autofix_permuted.log", help="Permuted-match sidecar log")
    ap.add_argument("--clear_log", action="store_true",
                    help="Truncate the permuted log after a successful fix "
                         "so entries are not applied twice on the next run")
    args = ap.parse_args()

    entries = parse_permuted_log(args.log)
    if not entries:
        print("[fix_permuted] No entries in %s — nothing to do." % args.log)
        return

    bak = args.pack + ".prefixbak"
    try:
        open(bak).close()
    except IOError:
        import shutil as _sh
        _sh.copy2(args.pack, bak)
        print("[fix_permuted] Backup created: %s" % bak)

    with open(args.pack, "r") as f:
        lines = f.readlines()

    fixed = 0
    for entry in entries:
        kind = entry["kind"]
        if kind not in FIXERS:
            # bonds are always symmetric — skip
            continue
        fixer = FIXERS[kind]
        fixer(lines, **entry)
        print("[fix_permuted] Fixed %s type %d ('%s' matched as '%s')"
              % (kind, entry["new_type"], entry["label"], entry["matched"]))
        fixed += 1

    with open(args.pack, "w") as f:
        f.writelines(lines)

    print("[fix_permuted] Done — %d type(s) corrected in %s." % (fixed, args.pack))

    if args.clear_log:
        open(args.log, "w").close()
        print("[fix_permuted] Cleared %s." % args.log)


if __name__ == "__main__":
    main()
