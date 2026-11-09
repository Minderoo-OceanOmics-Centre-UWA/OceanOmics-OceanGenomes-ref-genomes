#!/usr/bin/env python3
"""Compile what decontamination removed from each haplotype, for ref_decontamination.

One row per assembly and haplotype (hap1, hap2), from the DECONTAMINATION subworkflow's outputs:
  input_seqs, input_bp   the YaHS scaffolds FCS-GX was given (read from the scaffold FASTA)
  fcs_*                  FCS-GX report rows by action; bp = end - start + 1
  tiara_*                the contigs TIARA_TIARA removes (modules/nf-core/tiara/tiara/main.nf):
                         mitochondrion, plastid, and prokarya / bacteria / archaea as prokaryote,
                         with bp from the scaffold FASTA
Not the module's tiara_filter_summary.txt: its bp column is always 0.

Sources:
  scratch (default)  <base>/OG*/05-decontamination/{NCBI/out,tiara}/ and 04-scaffolding/yahs/
  --acacia REMOTE    the backups, REMOTE/OG*/OG*_v<date>.<version>/{decontamination/{NCBI,tiara},assembly}/,
                     read with `rclone cat`: the FASTAs are streamed, nothing is kept

Writes final_decontam_report.tsv for 06a_push_decontam_results_to_sqldb.py. Assemblies it cannot
read are listed on stderr and left out; a haplotype whose Tiara contigs are not all in the FASTA
gets no bp (and is listed).

    python3 06_decontam_compile.py [--base DIR | --acacia pawsey0964:oceanomics-refassemblies] [--only OG2089,OG90]
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

KEY = r"(?P<og>OG\d+)_v?(?P<date>\d{6})\.(?P<ver>(?:hic|hifi|ont)\d+)"
TIARA_RE = re.compile(rf"^{KEY}\.(?P<stage>\d)\.yahs\.(?P<hap>hap[12])\.tiara\.txt$")
FCS_RE = re.compile(rf"^{KEY}_?\.(?P<stage>\d)\.yahs\.(?P<hap>hap[12])\.NCBI\.\d+\.fcs_gx_report\.txt$")
FASTA_RE = re.compile(rf"^{KEY}\.(?P<stage>\d)\.yahs\.(?P<hap>hap[12])(?P<suffix>_scaffolds_final|_scaffolds)?\.fa(?:sta)?$")
FASTA_PREFERENCE = ("_scaffolds_final", "_scaffolds", None)

MITO, PLASTID, PROKARYOTE = "mito", "plastid", "prokaryote"
FCS_ACTIONS = ("EXCLUDE", "TRIM", "FIX", "REVIEW")
COLUMNS = ["og_id", "seq_date", "version", "haplotype", "input_seqs", "input_bp",
           "fcs_exclude_seqs", "fcs_exclude_bp", "fcs_trim_seqs", "fcs_trim_bp", "fcs_fix_seqs",
           "fcs_review_seqs", "tiara_mito_seqs", "tiara_mito_bp", "tiara_plastid_seqs",
           "tiara_plastid_bp", "tiara_prokaryote_seqs", "tiara_prokaryote_bp"]


# -- reading files --------------------------------------------------------------------------------

class Source:
    """Files of one place: a scratch directory, or an rclone remote."""

    def __init__(self, base=None, remote=None):
        self.base, self.remote = base, remote

    def list(self):
        """Relative paths of every candidate file."""
        if self.remote:
            out = subprocess.run(
                ["rclone", "lsf", "-R", "--fast-list", "--files-only",
                 "--include", "*/decontamination/**", "--include", "*/assembly/*yahs*", self.remote],
                check=True, capture_output=True, text=True).stdout
            return [p for p in out.splitlines() if p]
        paths = []
        for og in sorted(Path(self.base).glob("OG*")):
            for sub in ("05-decontamination/NCBI/out", "05-decontamination/tiara", "04-scaffolding/yahs"):
                d = og / sub
                if d.is_dir():
                    paths += [str(p.relative_to(self.base)) for p in sorted(d.iterdir()) if p.is_file()]
        return paths

    def lines(self, rel):
        if self.remote:
            proc = subprocess.Popen(["rclone", "cat", f"{self.remote}/{rel}"], stdout=subprocess.PIPE,
                                    text=True, bufsize=1 << 20)
            yield from proc.stdout
            if proc.wait() != 0:
                raise OSError(f"rclone cat {self.remote}/{rel} exited {proc.returncode}")
        else:
            with open(os.path.join(self.base, rel)) as f:
                yield from f


# -- parsing --------------------------------------------------------------------------------------

def tiara_removed(lines):
    """{contig: class} for the contigs the TIARA_TIARA module removes. Its rule is grep -w of each
    class word over the report line; the words can only appear in the two class columns."""
    removed = {}
    for i, line in enumerate(lines):
        cols = line.rstrip("\n").split("\t")
        if i == 0 or len(cols) < 3:
            continue
        classes = {cols[1], cols[2]}
        if "mitochondrion" in classes:
            removed[cols[0]] = MITO
        elif "plastid" in classes:
            removed[cols[0]] = PLASTID
        elif classes & {"prokarya", "bacteria", "archaea"}:
            removed[cols[0]] = PROKARYOTE
    return removed


def fcs_actions(lines):
    """{action: (distinct sequences, bp)} from an FCS-GX report."""
    seqs, bp = defaultdict(set), defaultdict(int)
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        seq_id, start, end, action = cols[0], int(cols[1]), int(cols[2]), cols[4]
        seqs[action].add(seq_id)
        bp[action] += end - start + 1
    return {a: (len(seqs[a]), bp[a]) for a in seqs}


def fasta_lengths(lines, wanted):
    """(sequence count, total bp, {id: length} for the ids in wanted)."""
    n, total, lengths = 0, 0, {}
    current, length = None, 0
    for line in lines:
        if line.startswith(">"):
            if current is not None:
                total += length
                if current in wanted:
                    lengths[current] = length
            n += 1
            current, length = line[1:].split()[0], 0
        else:
            length += len(line.strip())
    if current is not None:
        total += length
        if current in wanted:
            lengths[current] = length
    return n, total, lengths


def compile_row(key, tiara_lines, fcs_lines, fasta_lines):
    og, date, ver, hap = key
    removed = tiara_removed(tiara_lines)
    actions = fcs_actions(fcs_lines)
    row = {"og_id": og, "seq_date": date, "version": ver, "haplotype": hap}
    for action in FCS_ACTIONS:
        n, bp = actions.get(action, (0, 0))
        row[f"fcs_{action.lower()}_seqs"] = n
        if action in ("EXCLUDE", "TRIM"):
            row[f"fcs_{action.lower()}_bp"] = bp
    for cls in (MITO, PLASTID, PROKARYOTE):
        row[f"tiara_{cls}_seqs"] = sum(1 for c in removed.values() if c == cls)
    problem = None
    if fasta_lines is None:
        row["input_seqs"] = row["input_bp"] = None
        problem = "no scaffold FASTA"
    else:
        n, total, lengths = fasta_lengths(fasta_lines, set(removed))
        row["input_seqs"], row["input_bp"] = n, total
        missing = set(removed) - set(lengths)
        if missing:
            problem = f"{len(missing)} Tiara contig(s) not in the scaffold FASTA"
    for cls in (MITO, PLASTID, PROKARYOTE):
        row[f"tiara_{cls}_bp"] = (None if problem else
                                  sum(lengths[c] for c, k in removed.items() if k == cls))
    return row, problem


# -- matching files to assemblies -----------------------------------------------------------------

def index(paths):
    """{(og, date, version, hap): {'tiara': path, 'fcs': path, 'fasta': path}}; stage 2 files win
    over stage 1 ones, and _scaffolds_final over _scaffolds over a bare .fa."""
    found = defaultdict(dict)
    for rel in paths:
        name = rel.rsplit("/", 1)[-1]
        for kind, rx in (("tiara", TIARA_RE), ("fcs", FCS_RE), ("fasta", FASTA_RE)):
            m = rx.match(name)
            if not m:
                continue
            key = (m["og"], m["date"], m["ver"], m["hap"])
            rank = (-int(m["stage"]),) + ((FASTA_PREFERENCE.index(m["suffix"]),) if kind == "fasta" else ())
            if kind not in found[key] or rank < found[key][kind][0]:
                found[key][kind] = (rank, rel)
    return {k: {kind: v[1] for kind, v in kinds.items()} for k, kinds in found.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--base", default=f"/scratch/pawsey0964/{os.environ.get('USER', '')}/ref-gen")
    src.add_argument("--acacia", metavar="REMOTE", help="e.g. pawsey0964:oceanomics-refassemblies")
    ap.add_argument("--only", help="comma-separated OG IDs")
    ap.add_argument("--out", default="final_decontam_report.tsv")
    args = ap.parse_args()

    source = Source(remote=args.acacia) if args.acacia else Source(base=args.base)
    only = set(args.only.split(",")) if args.only else None
    files = index(source.list())
    rows, skipped = [], []
    for key in sorted(files):
        if only and key[0] not in only:
            continue
        f = files[key]
        if "tiara" not in f or "fcs" not in f:
            skipped.append((key, "no " + " and ".join(k for k in ("tiara", "fcs") if k not in f)))
            continue
        try:
            row, problem = compile_row(key, source.lines(f["tiara"]), source.lines(f["fcs"]),
                                       source.lines(f["fasta"]) if "fasta" in f else None)
        except (OSError, ValueError, IndexError) as exc:
            skipped.append((key, f"unreadable: {exc}"))
            continue
        rows.append(row)
        if problem:
            skipped.append((key, f"{problem}; bp left blank"))
        print(f"{'_'.join(key[:2])}.{key[2]} {key[3]}: in {row['input_seqs']}, "
              f"Tiara {row['tiara_mito_seqs']}/{row['tiara_plastid_seqs']}/{row['tiara_prokaryote_seqs']}, "
              f"FCS excl {row['fcs_exclude_seqs']} trim {row['fcs_trim_seqs']}", file=sys.stderr, flush=True)

    with open(args.out, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join("" if row[c] is None else str(row[c]) for c in COLUMNS) + "\n")
    print(f"\n{len(rows)} haplotype row(s) written to {args.out}", file=sys.stderr)
    for key, why in skipped:
        print(f"  not complete: {'_'.join(key[:2])}.{key[2]} {key[3]}: {why}", file=sys.stderr)


if __name__ == "__main__":
    main()
