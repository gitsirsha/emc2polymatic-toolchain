# emc2polymatic-toolchain
Enhanced Monte Carlo (EMC) and Polymatic workflow streamlining for Polymer Building

Contributors: <names..>

## Why this exists

In our group, we often need to build **amorphous polymer systems**, especially **glassy/high-Tg polymers** and **rigid monomers** (e.g., PIM-like chemistries). Two common tools in this space are **EMC** and **Polymatic**—both powerful, but with different strengths and practical limitations.

### EMC: fast monomer → typed topology → build (great for many cases)
EMC is extremely convenient when you want to go from **SMILES/chemistry + force field choice + density/box targets** to a **typed structure/topology** quickly. It is also valuable because it supports several force fields (e.g., PCFF/COMPASS/CHARMM/OPLS, etc.) and can generate simulation inputs/exports for common engines.

However, for some **multifunctional monomers** (e.g., “2=2” reactive connections needed for certain rigid/glassy monomers), the default polymerization/build assumptions can be limiting. In practice, EMC workflows are often most straightforward when the connectivity growth resembles **single-bond connections between monomers**.

### Polymatic: packing + flexible polymerization rules + strong equilibration
Polymatic is very useful for:
- **packing + polymerization** workflows where you need more flexible bonding schemes (including cases that can emulate “2=2” style connections with the right setup),
- and its widely used **multi-stage equilibration protocol** (commonly referenced as a ~21-step workflow) that tends to work well for **glassy/high-Tg systems**.

A common bottleneck is that Polymatic typically assumes you can provide (or derive) the **monomer topology/typing** in a form it can consume reliably—something EMC can often produce more conveniently for many chemistries.

### The idea: use each tool for what it’s best at
This repository aims to make the **EMC → Polymatic handoff** smoother and more reproducible:
- consistent folder structure and “recipes”
- minimal configs/templates
- clear, repeatable steps and checks
- a practical bridge for group members who want a reliable starting point without re-learning both ecosystems every time

## Project status

**Ongoing / under active development.**

This repo currently focuses on documentation + templates. Automation/CLI glue will be added iteratively as workflows stabilize.

Planned milestones:
- [ ] Minimal end-to-end example (“toy” system) documenting EMC → Polymatic → LAMMPS outputs
- [ ] Template library for common monomers / force fields used in the group
- [ ] Optional wrapper scripts (one-command runs) + validation checks
- [ ] Troubleshooting guide for common failure modes (typing/connectivity/packing)

## Citations

### EMC
In any publication of scientific results based in part or completely on the use of EMC, please include:

P.J. in ’t Veld and G.C. Rutledge, *Macromolecules* **2003**, 36, 7358.

### Polymatic
L.J. Abbott, K.E. Hart, and C.M. Colina, *Theoretical Chemistry Accounts* **2013**, 132:1334. DOI: 10.1007/s00214-013-1334-z

## Third-party notice
This repository does **not** include EMC or Polymatic. They must be installed separately and are governed by their own licenses/citation requirements.
This project is not affiliated with the original authors/maintainers of EMC or Polymatic.

