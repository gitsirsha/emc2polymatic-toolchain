#!/blue/jsampath/ganguly.sirsha/.conda/envs/polymer-building/bin/python <-- edit this path from python of your environment

################################################################################
#
# swap_improper_comments.py
#
# Author: Sirsha Ganguly
#
# Description:
#   Preprocesses an EMC-generated PCFF LAMMPS data file (pack.lmps / data.lmps)
#   to fix the improper-type comment convention before running Polymatic.
#
#   Although in the topology section impropers are written by EMC as:
# 
#     Impropers: a,CENTER,b,c (correct Class-2 convention also used by Polymatic)
#  
#   EMC writes force-field parameters with the center atom at index 0
#   in the type label comments (the one with #):
#
#     Improper Coeffs:   CENTER,a,b,c    (EMC / force-field convention)
#     AngleAngle Coeffs: CENTER,a,b,c
#
#   Polymatic (polym.pl + Polymatic.pm) and LAMMPS class2 both expect the
#   center atom at index 1 (second position) in the topology:
#
#     Impropers block:   a,CENTER,b,c    (LAMMPS topology convention)
#
#   autofix_from_log_improper.py writes new type comments in the topology
#   convention (center at index 1) to stay consistent with what polym.pl
#   reports in its error messages.  For the initial entries already present
#   in the file (written by EMC) to match the same convention, this script
#   must be run ONCE on the pack file before starting the polymerization loop.
#
# Usage:
#   python swap_improper_comments.py -i pack.lmps [-o pack_fixed.lmps]
#
#   If -o is omitted the input file is edited in place (a .bak backup is
#   created automatically).
#
# Sections modified:
#   - Improper Coeffs
#   - AngleAngle Coeffs
#
################################################################################

import argparse
import re
import shutil


SWAP_SECTIONS = ["Improper Coeffs", "AngleAngle Coeffs"]


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


def swap_comment_center(comment):
    """
    Given a comma-separated 4-token label, move index 0 to index 1:
      'CENTER,a,b,c'  ->  'a,CENTER,b,c'
    Returns the original string unchanged if it is not a 4-token label.
    """
    toks = [t.strip() for t in comment.split(",")]
    if len(toks) != 4:
        return comment
    center, a, b, c = toks
    return ",".join([a, center, b, c])


def process_lines(lines):
    changed = 0
    for sec_header in SWAP_SECTIONS:
        sec = find_section(lines, sec_header)
        if sec is None:
            continue
        _, b0, b1 = sec
        for idx in range(b0, b1):
            ln = lines[idx]
            if "#" not in ln:
                continue
            left, right = ln.split("#", 1)
            old_comment = right.strip()
            new_comment = swap_comment_center(old_comment)
            if new_comment != old_comment:
                lines[idx] = left + "# " + new_comment + "\n"
                changed += 1
    return changed


def main():
    ap = argparse.ArgumentParser(
        description="Swap improper-type comment convention: center-first -> center-second."
    )
    ap.add_argument("-i", "--input",  required=True,  help="Input LAMMPS data file")
    ap.add_argument("-o", "--output", default=None,   help="Output file (default: edit in place)")
    args = ap.parse_args()

    with open(args.input, "r") as f:
        lines = f.readlines()

    if args.output is None:
        bak = args.input + ".bak"
        shutil.copy2(args.input, bak)
        print("[swap] Backup: %s" % bak)
        out_path = args.input
    else:
        out_path = args.output

    changed = process_lines(lines)

    with open(out_path, "w") as f:
        f.writelines(lines)

    print("[swap] %d comment(s) updated in %s" % (changed, out_path))


if __name__ == "__main__":
    main()
