#Author, date: Laura, Ozan, and Claude, Sept 1 2026
#Motivation: Process the OGSharedness table
#Intent: Filters the table and converts # seqs to pres/abs and allows users to specify targets
#Dependencies: pandas
#Inputs: OG sharedness spreadsheet  
#Outputs: Filtered OGSharedness table or Presence/absence spreadsheets
#Example: usage to run are below

"""
Filter the OGSharedness CSV or convert into a presence/absence CSV (Species, MinorClades, 
MajorClades and the major-clade sums recomputed from presence). Sequences and Paralogness 
are copied from the input. Species codes are XX_yy_Zzzz (major = XX, minor = XX_yy).

Two output modes:
    default         presence/absence output, as described above.
    --filter_only   filter the original file only: no presence/absence
                    conversion. Output = the original leading columns
                    (everything before the first XX_yy_Zzzz column), copied
                    unchanged, plus the filtered species columns with their
                    original sequence counts. The filters below, plus
                    --target_clade and --non_target, still apply.

Arguments:
  -i, --input        Input OGSharedness.csv (required)
  -o, --output       Output CSV name (default: auto-named from the filters and thresholds used)

  --target_clade     Keep only OGs that meet a condition on a 2, 5 or 10-character code,
                     e.g. 'Am_tu>5' (more than 5 Am_tu species present).
                     Operators: > >= < <= == !=
                     Several conditions (separated by space or comma) must ALL be true (AND).
                     Use '|' inside one condition for OR, e.g. 'Am>0|Ex_he>0'.
                     Counts always use ALL species, even with --filter_taxa or --non_target.

  --non_target       2, 5 or 10-character code(s) for species columns to drop from the output.
                     These species are still counted in all summaries and --target_clade conditions.

  --filter_taxa      Show only species columns starting with these 2, 5 or 10-character codes,
                     e.g. Am, Am_tu, Am_tu_Aint. Leading columns are always kept.

  --filter_OG        Keep only these OGs. Give a file (one OG ID per line, '#' lines skipped)
                     or type OG IDs directly.

  --filter_only      Keep the original sequence counts instead of converting to presence/absence.
                     Output filename starts with 'Filtered' instead of 'PresAbs'.

Example:
  python script.py -i OGSharedness.csv --target_clade 'Am>0|Ex_he>0' Op>=3 --non_target Ba,Za


Output columns (presence/absence mode):
    OG, Sequences, Species, Paralogness, MinorClades, MajorClades,
    Am, EE, Ex, Op, Pl, Sr, Euks, <blank>, Ba, Za, Proks, <blank>,
    [one count column per non-major --target_clade], <blank>,
    <species presence/absence columns>

Output filename (unless -o is given) encodes run specific values and
filters, so different runs never overwrite each other:
    PresAbs.csv
    PresAbs_Am_tu-gt5.csv                 (--target_clade "Am_tu>5")
    PresAbs_Am_tu-gt5_Ba-ge3.csv          (--target_clade "Am_tu>5" "Ba>=3")
    PresAbs_Am_tu-gt5_no-Op-Pl.csv        (... --non_target Op Pl)
    PresAbs_Am-gt0-or-Ex_he-gt0.csv       (--target_clade "Am>0|Ex_he>0")
    PresAbs_taxa-Am-Op_OGs-mylist.csv     (--filter_taxa Am Op --filter_OG mylist.txt)
    Filtered_taxa-Am-Op_OGs-mylist.csv    (same, with --filter_only)
  (OG IDs typed directly are named individually if <= 3, else as "<N>OGs")
Species absent from every kept OG are dropped and listed in
<output>_removed_species.txt.

Usage:
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv -o custom_name.csv
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv --target_clade "Am_tu>5" "Ba>3"
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv --target_clade "Am_tu>5" --non_target Op Pl
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv --target_clade "Am>0|Ex_he>0" "Ba>2"
        (--target_clade logic: separate conditions, by space or comma, are ANDed;
         alternatives joined by '|' inside one condition are ORed. The example
         keeps OGs with (Am>0 OR Ex_he>0) AND Ba>2.)
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv -filter_taxa Am Op_me Pl_gr_Atha
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv -filter_OG my_OGs.txt --target_clade "Am>3"
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv -filter_OG OG6_129803 OG6_170434
    python Target_Pres_Abs_from_OGSharedness.py -i input.csv -filter_only -filter_taxa Am -filter_OG my_OGs.txt
"""

import argparse
import operator
import os
import re
import sys

import pandas as pd

META_COLS = ["OG", "Sequences", "Species", "Paralogness", "MinorClades", "MajorClades"]
EUK_MAJORS = ["Am", "EE", "Ex", "Op", "Pl", "Sr"]
PROK_MAJORS = ["Ba", "Za"]
MAJOR_CLADE_COLS = EUK_MAJORS + PROK_MAJORS
SPECIES_RE = re.compile(r"^[A-Za-z0-9]{2}_[A-Za-z0-9]{2}_[A-Za-z0-9]{4}$")
OPS = {">=": (operator.ge, "ge"), "<=": (operator.le, "le"), "==": (operator.eq, "eq"),
       "!=": (operator.ne, "ne"), ">": (operator.gt, "gt"), "<": (operator.lt, "lt")}
TARGET_RE = re.compile(r"^\s*(\S+?)\s*(>=|<=|==|!=|>|<)\s*(-?\d+(?:\.\d+)?)\s*$")
    

# Function for parsing multiple target-clade inputs (OR statement)
def parse_target_group(expr):
    """'Am>0|Ex_he>0' -> list of OR alternatives, plus a filename tag."""
    alts = [parse_target(a) for a in expr.split("|") if a.strip()]
    if not alts:
        sys.exit(f"Empty --target_clade condition: {expr!r}")
    return alts, "-or-".join(a[3] for a in alts)


# Function for parsing target-clade input
def parse_target(expr):
    m = TARGET_RE.match(expr)
    if not m:
        sys.exit(f"Could not parse --target_clade {expr!r}; expected e.g. 'Am_tu>5'.")
    clade, sym, raw = m.groups()
    func, tag = OPS[sym]
    return clade, func, float(raw), f"{clade}-{tag}{raw}"

# Function for turning comma separated argument values into a flat list (e.g., --filter_taxa "Am ,Am_tu")
def flatten(groups):
    """Flatten argparse's list-of-lists, split commas, drop duplicates (order kept)."""
    items = (p.strip() for g in groups or [] for item in g for p in item.split(","))
    return list(dict.fromkeys(p for p in items if p))

# Function for reading OG list
def read_og_list(values):
    """Each value is a file (one OG per line, first column, '#' lines skipped)
    or an OG ID. Returns (unique OG IDs, filename tags)."""
    ogs, tags, inline = [], [], []
    for v in flatten(values):
        if os.path.isfile(v):
            with open(v) as fh:
                ogs += [re.split(r"[\s,]+", l.strip())[0] for l in fh
                        if l.strip() and not l.startswith("#")]
            tags.append(os.path.splitext(os.path.basename(v))[0])
        else:
            inline.append(v)
    ogs += inline
    if inline:
        tags += inline if len(inline) <= 3 else [f"{len(inline)}OGs"]
    return list(dict.fromkeys(ogs)), tags

# Function for filtering columns
def filter_taxa_columns(species, codes):
    """Keep species columns starting with any 2-, 5- or 10-character code."""
    bad = [c for c in codes if len(c) not in (2, 5, 10)]
    if bad:
        sys.exit(f"--filter_taxa codes must be 2, 5 or 10 characters (e.g. Am, Am_tu, Am_tu_Aint): {bad}")
    for c in codes:
        if not any(s.startswith(c) for s in species):
            print(f"Warning: filter_taxa '{c}' did not match any species columns.", file=sys.stderr)
    kept = [s for s in species if s.startswith(tuple(codes))]
    if not kept:
        sys.exit("No species columns left after --filter_taxa.")
    return kept

# Function for writing the output file
def build_output_filename(tags, non_target, taxa=(), og_tags=(), prefix="PresAbs"):
    name = prefix + "".join(f"_{t}" for t in tags)
    if non_target:
        name += "_no-" + "-".join(non_target)
    if taxa:
        name += "_taxa-" + "-".join(taxa)
    if og_tags:
        name += "_OGs-" + "-".join(og_tags)
    return name + ".csv"

# Main function
def main():
    # Parse arguments 
    ap = argparse.ArgumentParser(description="Filters an OGSharedness table or converts it into a presence/absence CSV")
    ap.add_argument("-i", "--input", dest="input_csv", required=True, help="Input OGSharedness.csv")
    ap.add_argument("-o", "--output", help="Output CSV (default: auto-named based on filters and threshold values)")
    ap.add_argument("--target_clade", nargs="+", action="append",
                    help="2, 5, or 10 digit codes followed by condition, e.g. 'Am_tu>5'. "
                         "Conditions: > >= < <= == !=. Separate conditions (space or comma) are processed as AND statement; "
                         "use '|' inside one condition for OR, e.g. 'Am>0|Ex_he>0'.")
    ap.add_argument("--non_target", nargs="+", action="append",
                    help="2,5, or 10 digit code for species to be filtered out (still counted in all summaries)")
    ap.add_argument("-filter_taxa", "--filter_taxa", nargs="+", action="append",
                    help="Keep only species columns starting with these 2, 5 or 10-digit codes, "
                         "e.g. Am (all Am species), Am_tu, Am_tu_Aint. Leading columns are always kept. "
                         "All counts and --target_clade still use ALL species, so they match the "
                         "original file.")
    ap.add_argument("-filter_OG", "--filter_OG", nargs="+", action="append",
                    help="Keep only these OGs. Each value is a file with one OG ID per line "
                         "(first column used, '#' lines skipped) or an OG ID typed directly.")
    ap.add_argument("-filter_only", "--filter_only", action="store_true",
                    help="Filter the original file only, keeping the original sequence counts (no "
                         "presence/absence conversion). Output contains original leading columns (copied "
                         "unchanged) and filtered species columns. Filename starts with "
                         "'Filtered' instead of 'PresAbs'.")
    args = ap.parse_args()

    df = pd.read_csv(args.input_csv)
    missing = [c for c in META_COLS if c not in df.columns]
    if missing:
        sys.exit(f"Input file is missing expected columns: {missing}")
    known = set(META_COLS) | set(MAJOR_CLADE_COLS)
    species = [c for c in df.columns if c not in known and SPECIES_RE.match(c)]
    ignored = [c for c in df.columns if c not in known and c not in species]
    if ignored:
        print(f"Warning: ignoring columns that are not species codes: {ignored}", file=sys.stderr)
    if not species:
        sys.exit("No species columns detected in input file.")
    all_species = species  # every count/summary/target uses ALL species

    # Optional filters: -filter_taxa only chooses which species COLUMNS are shown;
    # -filter_OG chooses rows. Neither changes any per-OG count.
    filter_taxa = flatten(args.filter_taxa)
    if filter_taxa:
        species = filter_taxa_columns(species, filter_taxa)
        print(f"--filter_taxa: keeping {len(species)} species columns.")
    og_list, og_tags = read_og_list(args.filter_OG)
    if og_list:
        not_found = set(og_list) - set(df["OG"].astype(str))
        if not_found:
            print(f"Warning: {len(not_found)} OG(s) from --filter_OG not in input, e.g. "
                  f"{sorted(not_found)[:5]}", file=sys.stderr)
        df = df[df["OG"].astype(str).isin(og_list)]
        if df.empty:
            sys.exit("No OGs left after --filter_OG.")
        print(f"--filter_OG: keeping {len(df)} of {len(og_list)} listed OGs.")

    groups = [parse_target_group(e) for e in flatten(args.target_clade)]
    targets = [alt for alts, _ in groups for alt in alts]  # every single condition
    non_target = tuple(flatten(args.non_target))
    for p in non_target:
        if not any(c.startswith(p) for c in species):
            print(f"Warning: non_target '{p}' did not match any species columns.", file=sys.stderr)

    # Presence/absence and clade counts, always computed from ALL species
    pres = (df[all_species].apply(pd.to_numeric, errors="coerce").fillna(0) > 0).astype(int)
    by_minor = pres.T.groupby(lambda c: c[:5]).sum().T
    by_major_all = pres.T.groupby(lambda c: c[:2]).sum().T
    major = by_major_all.reindex(columns=MAJOR_CLADE_COLS, fill_value=0)

    def count(prefix):
        cols = [c for c in all_species if c.startswith(prefix)]
        if not cols:
            sys.exit(f"target_clade '{prefix}' matches no species columns (codes are case-sensitive).")
        return pres[cols].sum(axis=1)

    mask = pd.Series(True, index=df.index)
    for alts, _ in groups:  # AND across conditions, OR within a '|' condition
        any_alt = pd.Series(False, index=df.index)
        for clade, func, value, _ in alts:
            any_alt |= func(count(clade), value)
        mask &= any_alt

    # Species columns: drop non-targets, then species absent from every kept OG
    candidates = [c for c in species if not c.startswith(non_target)]
    present = pres.loc[mask, candidates].any()
    removed = list(present.index[~present])
    kept = list(present.index[present])
    if non_target:
        print(f"Dropped {len(species) - len(candidates)} non-target species columns (of the kept taxa).")

    if args.filter_only:
        # Keep original leading columns (before the first species code) + original counts
        first = next(i for i, c in enumerate(df.columns) if SPECIES_RE.match(c))
        raw = pd.read_csv(args.input_csv, dtype=str, keep_default_na=False)  # exact original text
        out = raw.loc[mask.index[mask], list(df.columns[:first]) + kept]
    else:
        blank = pd.Series("", index=df.index)
        cols = [("OG", df["OG"]), ("Sequences", df["Sequences"]), ("Species", pres.sum(axis=1)),
                ("Paralogness", df["Paralogness"]),
                ("MinorClades", (by_minor > 0).sum(axis=1)), ("MajorClades", (by_major_all > 0).sum(axis=1))]
        cols += [(m, major[m]) for m in EUK_MAJORS] + [("Euks", major[EUK_MAJORS].sum(axis=1)), ("", blank)]
        cols += [(m, major[m]) for m in PROK_MAJORS] + [("Proks", major[PROK_MAJORS].sum(axis=1)), ("", blank)]
        extra = [c for c in dict.fromkeys(t[0] for t in targets) if c not in MAJOR_CLADE_COLS]
        if extra:  # suffix full species codes so they don't clash with that species' own column
            cols += [(f"{c}_count" if c in all_species else c, count(c)) for c in extra] + [("", blank)]
        cols += [(c, pres[c]) for c in kept]
        out = pd.concat([s for _, s in cols], axis=1).loc[mask]
        out.columns = [n for n, _ in cols]

    output = args.output or build_output_filename([tag for _, tag in groups], non_target, filter_taxa, og_tags,
                                                  "Filtered" if args.filter_only else "PresAbs")
    if os.path.abspath(output) == os.path.abspath(args.input_csv):
        sys.exit("Output path is the same as the input file; exiting to not overwrite.")
    if os.path.exists(output):
        print(f"Warning: overwriting existing {output}", file=sys.stderr)
    out.to_csv(output, index=False)
    print(f"Wrote {len(out)} rows to {output}")

    if removed:
        print(f"Removed {len(removed)} all-absent species columns.")

if __name__ == "__main__":
    main()
