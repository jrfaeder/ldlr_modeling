"""
LDLR structural features from PDB 9BDE and PDB 1N7D
===================================================

PURPOSE
  Computes three per-residue structural features for LDLR and adds them to
  the variant table (each variant inherits the values of its position):
    dist_apob100   minimum Cα–Cα distance (Å) from the LDLR residue to any
                   ApoB100 residue                       [9BDE only]
    is_in_helix    1 if DSSP assigns a helix (H, G or I), else 0
    is_in_strand   1 if DSSP assigns a β-strand (E), else 0
    ss_source      structure the helix/strand values come from
                   ("9BDE" or "1N7D"); bookkeeping, not a model feature

STRUCTURES
  9BDE  LDLR bound to ApoB100, cryo-EM (Reimund et al., Nature 2025).
        LDLR = chain R, residues 66-354 (LA2 to EGF-A); ApoB100 = chain A.
        The other chains (B: maltose-binding protein fusion; H, L: antibody
        Fab; N: nanobody) are not LDLR and are not used.
  1N7D  LDLR extracellular domain at endosomal pH, X-ray (Rudenko et al.,
        Science 2002). LDLR = chain A, residues 65-714 (LA2 to EGF-C,
        including the full β-propeller). Contains two engineered
        glycosylation-site substitutions, N515Q and N657Q.

  How the two are combined:
    dist_apob100            9BDE only (1N7D contains no ApoB100).
    is_in_helix/strand      9BDE where the residue is resolved there (the
                            ligand-bound complex); otherwise 1N7D.
  In the 265 residues resolved in both structures, DSSP agrees for 92.5%
  of helix and 96.6% of strand assignments.
  Residues in neither structure (signal, LA1, O-sugar, TM, cytoplasmic
  tail) are left empty (NaN).

  LDLR numbering in both structures matches the Tabet et al. data
  (precursor numbering). Before writing anything, the script checks every
  resolved residue against ref_aa3 in the variant table (allowing only the
  two engineered 1N7D substitutions).

INPUT
  --pdb_9bde, --pdb_1n7d   structure files (mmCIF or PDB) from rcsb.org.
                           The same coordinates are in Tabet Data S8 (9BDE)
                           and Data S7 (1N7D), as PyMOL sessions.
  --csv                    variant table (needs position, ref_aa3)

OUTPUT
  --out           per-position table: position, dist_apob100, is_in_helix,
                  is_in_strand, ss_source, dssp
  --out_variants  (optional) the variant table with dist_apob100,
                  is_in_helix, is_in_strand and ss_source added/replaced

USAGE
  python ldlr_structure_features.py --pdb_9bde 9BDE.cif --pdb_1n7d 1N7D.cif \\
      --csv ldlr_variant_scores_v5.csv --out structure_by_position.csv \\
      --out_variants ldlr_variant_scores_v5.csv

REQUIREMENTS
  numpy, pandas, biopython, DSSP (mkdssp; e.g. `conda install -c salilab
  dssp` or `apt install dssp`)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import DSSP, MMCIFParser, PDBParser

HELIX_CODES = {"H", "G", "I"}
STRAND_CODES = {"E"}
ENGINEERED_1N7D = {515: ("ASN", "GLN"), 657: ("ASN", "GLN")}


def load_model(path: str):
    path = Path(path)
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} \
        else PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))[0]


def protein_residues(chain):
    """Standard amino-acid residues with a Cα atom (ions and waters dropped)."""
    return [r for r in chain if r.id[0] == " " and "CA" in r]


def check_identity(residues, ref: pd.Series, name: str, allowed=None):
    allowed = allowed or {}
    bad = []
    for r in residues:
        pos, resn = r.id[1], r.get_resname()
        if pos in ref.index and ref[pos] != resn and allowed.get(pos) != (ref[pos], resn):
            bad.append((pos, resn, ref[pos]))
    if bad:
        raise ValueError(f"{name}: residues do not match LDLR numbering: {bad[:10]}")
    print(f"{name}: residue identity check passed "
          f"({sum(r.id[1] in ref.index for r in residues)} positions)")


def dssp_codes(model, path: str, chain_id: str) -> dict:
    dssp = DSSP(model, path, dssp="mkdssp")
    codes = {res_id[1]: values[2]
             for (ch, res_id), values in dssp.property_dict.items()
             if ch == chain_id and res_id[0] == " "}
    if not any(c in HELIX_CODES | STRAND_CODES for c in codes.values()):
        raise RuntimeError(f"DSSP assigned no helix or strand in {path}. Use the "
                           "original RCSB file (exports may lack the header DSSP needs).")
    return codes


def compute(pdb_9bde: str, pdb_1n7d: str, ref: pd.Series) -> pd.DataFrame:
    # 9BDE: LDLR chain R, ApoB100 chain A
    m9 = load_model(pdb_9bde)
    ldlr9 = protein_residues(m9["R"])
    apob = protein_residues(m9["A"])
    check_identity(ldlr9, ref, "9BDE chain R")
    apob_ca = np.array([r["CA"].coord for r in apob])
    dist = {r.id[1]: float(np.min(np.linalg.norm(apob_ca - r["CA"].coord, axis=1)))
            for r in ldlr9}
    ss9 = dssp_codes(m9, pdb_9bde, "R")

    # 1N7D: LDLR chain A
    m1 = load_model(pdb_1n7d)
    ldlr1 = protein_residues(m1["A"])
    check_identity(ldlr1, ref, "1N7D chain A", allowed=ENGINEERED_1N7D)
    ss1 = dssp_codes(m1, pdb_1n7d, "A")

    pos9 = {r.id[1] for r in ldlr9}
    pos1 = {r.id[1] for r in ldlr1}
    rows = []
    for pos in sorted(pos9 | pos1):
        if pos in pos9:
            code, source = ss9.get(pos, "-"), "9BDE"
        else:
            code, source = ss1.get(pos, "-"), "1N7D"
        rows.append(dict(position=pos,
                         dist_apob100=round(dist[pos], 3) if pos in dist else np.nan,
                         is_in_helix=int(code in HELIX_CODES),
                         is_in_strand=int(code in STRAND_CODES),
                         ss_source=source, dssp=code))
    table = pd.DataFrame(rows)
    print(f"9BDE: {len(pos9)} LDLR residues; 1N7D: {len(pos1)}; "
          f"combined: {len(table)} positions "
          f"({(table.ss_source == '1N7D').sum()} from 1N7D only)")
    return table


def main():
    ap = argparse.ArgumentParser(description="LDLR structural features (9BDE + 1N7D)")
    ap.add_argument("--pdb_9bde", required=True)
    ap.add_argument("--pdb_1n7d", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", default="structure_by_position.csv")
    ap.add_argument("--out_variants", default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.csv, low_memory=False)
    ref = (df.dropna(subset=["ref_aa3"]).drop_duplicates("position")
             .set_index("position")["ref_aa3"].str.upper())

    table = compute(args.pdb_9bde, args.pdb_1n7d, ref)
    table.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")

    if args.out_variants:
        t = table.set_index("position")
        for col in ["dist_apob100", "is_in_helix", "is_in_strand", "ss_source"]:
            df[col] = df["position"].map(t[col])
        df.to_csv(args.out_variants, index=False)
        print(f"Wrote {args.out_variants}: dist_apob100 for "
              f"{df['dist_apob100'].notna().sum()} variants, secondary structure for "
              f"{df['is_in_helix'].notna().sum()} of {len(df)}")


if __name__ == "__main__":
    main()
