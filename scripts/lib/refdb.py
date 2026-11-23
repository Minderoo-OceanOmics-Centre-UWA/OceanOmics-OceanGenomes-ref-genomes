"""refdb: database writes for the reference-genome push scripts (Phase RG, schema v4).

The same file is vendored in OceanOmics-OceanGenomes-ref-genomes and
OceanOmics-OceanGenomes-post-curation as scripts/lib/refdb.py. Change both together:
tests/test_refdb_copies.py fails when they differ.

The tables are the relational design in genomes_db_review/docs/ref_genomes_runbook_rg0_rg5.md:
one row per assembly (ref_assembly), per part of it (ref_assembly_part: stage x haplotype), and
one table per tool keyed on the part. A tool that did not run has no row, not a row of NULLs.

Table names are unqualified, so search_path decides where they resolve: v4 during the Phase RG
window (PGOPTIONS='-c search_path=v4,public'), public after the RG.6 swap.

Every write runs in one transaction per input file, and any error exits non-zero, so a caller
such as post_pipeline.sh or post-hichifi-assembly.py sees the failure.
"""

import configparser
import math
import re
import sys
from pathlib import Path

# Columns each tool table accepts. Anything else is a programming error.
TOOL_TABLES = {
    "ref_gfastats": ("num_contigs", "contig_n50", "contig_n50_size_mb", "num_scaffolds", "scaffold_n50",
                     "scaffold_n50_size_mb", "largest_scaffold", "largest_scaffold_size_mb",
                     "total_scaffold_length", "total_scaffold_length_size_mb", "gc_content_percent",
                     "num_gaps"),
    "ref_busco": ("dataset", "complete", "single_copy", "multi_copy", "fragmented", "missing", "n_markers",
                  "internal_stop_codon_percent", "scaffold_n50_bus", "contigs_n50_bus", "percent_gaps",
                  "number_of_scaffolds"),
    "ref_merqury_qv": ("unique_k_mers_assembly", "k_mers_total", "qv", "error"),
    "ref_merqury_completeness": ("solid_k_mers", "total_k_mers", "completeness"),
    "ref_omnic": ("total", "total_unmapped", "total_single_sided_mapped", "total_mapped", "total_dups",
                  "total_nodups", "cis", "trans"),
    "ref_seqkit_stats": ("format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len"),
    "ref_chromosome_stats": ("num_chromosomes", "pct_assigned", "pct_no_super", "num_seq_no_super",
                             "max_len_no_super"),
    # genomes_db_review migration 083: written by 06a_push_decontam_results_to_sqldb.py (stage 2).
    "ref_decontamination": ("input_seqs", "input_bp", "fcs_exclude_seqs", "fcs_exclude_bp", "fcs_trim_seqs",
                            "fcs_trim_bp", "fcs_fix_seqs", "fcs_review_seqs", "tiara_mito_seqs",
                            "tiara_mito_bp", "tiara_plastid_seqs", "tiara_plastid_bp",
                            "tiara_prokaryote_seqs", "tiara_prokaryote_bp"),
}
READ_SET_TABLES = {
    "ref_genomescope": ("homozygosity", "heterozygosity", "genomesize", "repeatsize", "uniquesize",
                        "modelfit", "errorrate"),
    "ref_hifi_adapter_filter": ("contam_reads",),
}
HIFI_READ_QC_KEY = ("og_id", "tissue", "ext_type", "lib_code", "run_id")
HIFI_READ_QC_COLS = ("barcode", "barcode_quality", "hifi_reads", "hifi_read_length", "hifi_read_quality",
                     "hifi_yield", "polymerase_read_length")
HAPLOTYPES = ("hap1", "hap2", "dual", "p_ctg", "a_ctg")
SEQ_DATE_RE = re.compile(r"^[0-9]{6}$")
VERSION_RE = re.compile(r"^(hic|hifi|ont)[0-9]+$")


# -- input parsing --------------------------------------------------------------------------------

def load_db_config(config_file):
    """psycopg2.connect() keyword arguments from the [postgres] section of an .ini file."""
    path = Path(config_file).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Config file '{path}' does not exist.")
    config = configparser.ConfigParser()
    config.read(path)
    if not config.has_section("postgres"):
        raise ValueError("Missing [postgres] section in config file.")
    missing = [k for k in ("dbname", "user", "password", "host", "port") if not config.has_option("postgres", k)]
    if missing:
        raise ValueError(f"Missing {missing} in [postgres] section of config file.")
    return {
        "dbname": config.get("postgres", "dbname"),
        "user": config.get("postgres", "user"),
        "password": config.get("postgres", "password"),
        "host": config.get("postgres", "host"),
        "port": config.getint("postgres", "port"),
    }


def _blank(val):
    if val is None:
        return True
    if isinstance(val, float) and math.isnan(val):
        return True
    s = str(val).strip()
    return s == "" or s.lower() == "nan"


def parse_int(val):
    """'1,234' / 1234.0 / '1234' -> 1234; blank, NaN or unparseable -> None."""
    if _blank(val):
        return None
    try:
        return int(float(str(val).strip().replace(",", "")))
    except ValueError:
        return None


def parse_float(val):
    """'12.5%' / '1,234.5' -> float; blank, NaN or unparseable -> None."""
    if _blank(val):
        return None
    try:
        return float(str(val).strip().replace(",", "").replace("%", ""))
    except ValueError:
        return None


def clean(val):
    """A value psycopg2 can bind: numpy scalars to Python, NaN and pandas NA to None."""
    if val is None:
        return None
    if hasattr(val, "item") and not isinstance(val, (str, bytes)):
        val = val.item()
    if isinstance(val, float) and math.isnan(val):
        return None
    if type(val).__name__ in ("NAType", "NaTType"):
        return None
    return val


def load_samplesheet(path, column):
    """{sample: value} for one samplesheet column, skipping blanks. None path -> {}."""
    import csv
    if not path:
        return {}
    out = {}
    with open(Path(path).expanduser(), newline="") as f:
        for row in csv.DictReader(f):
            og = (row.get("sample") or "").strip()
            val = (row.get(column) or "").strip()
            if og and val:
                out[og] = val
    return out


# -- assemblies and parts -------------------------------------------------------------------------

def check_part_key(og_id, seq_date, version, stage, haplotype):
    """Raise ValueError naming the bad part of a key, before the database sees it."""
    problems = []
    if not og_id or not str(og_id).startswith("OG"):
        problems.append(f"og_id {og_id!r}")
    if seq_date is None or not SEQ_DATE_RE.match(str(seq_date)):
        problems.append(f"seq_date {seq_date!r} (want YYMMDD)")
    if version is None or not VERSION_RE.match(str(version)):
        problems.append(f"version {version!r} (want hicN, hifiN or ontN)")
    if stage is None or str(stage) not in ("0", "1", "2", "3"):
        problems.append(f"stage {stage!r} (want 0-3)")
    if haplotype not in HAPLOTYPES:
        problems.append(f"haplotype {haplotype!r} (want one of {', '.join(HAPLOTYPES)})")
    if problems:
        raise ValueError(f"bad assembly key for {og_id}: " + "; ".join(problems))


def ensure_assembly(cur, og_id, seq_date, version):
    """assembly_id for the key, created if missing. ref_assembly.read_set_id (NOT NULL, migration 093)
    is filled by the ref_assembly_read_set_default trigger: a re-push keeps the assembly's read set,
    a new assembly gets the one read_set_runs() would give its date (created if missing)."""
    cur.execute(
        "INSERT INTO ref_assembly (og_id, seq_date, version) VALUES (%s, %s, %s) "
        "ON CONFLICT (og_id, seq_date, version) DO UPDATE SET og_id = EXCLUDED.og_id "
        "RETURNING assembly_id",
        (og_id, str(seq_date), version))
    return cur.fetchone()[0]


def ensure_part(cur, assembly_id, stage, haplotype):
    cur.execute(
        "INSERT INTO ref_assembly_part (assembly_id, stage, haplotype) VALUES (%s, %s, %s) "
        "ON CONFLICT (assembly_id, stage, haplotype) DO UPDATE SET stage = EXCLUDED.stage "
        "RETURNING part_id",
        (assembly_id, int(stage), haplotype))
    return cur.fetchone()[0]


def find_part(cur, og_id, seq_date, version, stage, haplotype):
    cur.execute(
        "SELECT p.part_id, a.assembly_id FROM ref_assembly a "
        "JOIN ref_assembly_part p ON p.assembly_id = a.assembly_id "
        "WHERE a.og_id = %s AND a.seq_date = %s AND a.version = %s AND p.stage = %s AND p.haplotype = %s",
        (og_id, str(seq_date), version, int(stage), haplotype))
    return cur.fetchone()


def prune(cur, part_id, assembly_id):
    """Drop a part with no tool rows left, then an assembly with no parts left: no empty rows."""
    tool_exists = " OR ".join(f"EXISTS (SELECT 1 FROM {t} x WHERE x.part_id = %(p)s)" for t in TOOL_TABLES)
    cur.execute(f"DELETE FROM ref_assembly_part WHERE part_id = %(p)s AND NOT ({tool_exists})", {"p": part_id})
    cur.execute("DELETE FROM ref_assembly a WHERE a.assembly_id = %s "
                "AND NOT EXISTS (SELECT 1 FROM ref_assembly_part p WHERE p.assembly_id = a.assembly_id)",
                (assembly_id,))


def upsert_metrics(cur, table, key, values):
    """Write one tool's values for one part.

    key: dict with og_id, seq_date, version, stage, haplotype.
    values: dict of that tool's columns. Columns not given are left as they are.
    If every value is None the tool's row for that part is removed (as the old wide-table upsert
    blanked the block), and an emptied part or assembly is removed with it.
    Returns 'upserted', 'deleted' or 'unchanged'.
    """
    if table not in TOOL_TABLES:
        raise ValueError(f"{table} is not a tool table")
    unknown = set(values) - set(TOOL_TABLES[table])
    if unknown:
        raise ValueError(f"{table} has no column(s) {sorted(unknown)}")
    k = {c: clean(key.get(c)) for c in ("og_id", "seq_date", "version", "stage", "haplotype")}
    check_part_key(**k)
    vals = {c: clean(v) for c, v in values.items()}

    if all(v is None for v in vals.values()):
        found = find_part(cur, **k)
        if not found:
            return "unchanged"
        part_id, assembly_id = found
        cur.execute(f"DELETE FROM {table} WHERE part_id = %s", (part_id,))
        prune(cur, part_id, assembly_id)
        return "deleted"

    assembly_id = ensure_assembly(cur, k["og_id"], k["seq_date"], k["version"])
    part_id = ensure_part(cur, assembly_id, k["stage"], k["haplotype"])
    cols = list(vals)
    quoted = [f'"{c}"' if c == "type" else c for c in cols]
    cur.execute(
        f"INSERT INTO {table} (part_id, {', '.join(quoted)}) "
        f"VALUES (%(part_id)s, {', '.join('%(' + c + ')s' for c in cols)}) "
        f"ON CONFLICT (part_id) DO UPDATE SET {', '.join(f'{q} = EXCLUDED.{q}' for q in quoted)}, loaded_at = now()",
        dict(vals, part_id=part_id))
    return "upserted"


# -- read sets and read QC ------------------------------------------------------------------------

# A samplesheet date can be the PacBio movie date, up to ~10 days off the run's seq_date (OG16:
# v240417, run PACB_240427_LA). Migration 062 dated the backfilled read sets with the same rule.
READ_SET_WINDOW_DAYS = 21


def read_set_runs(cur, og_id, sheet_date):
    """The OG's PacBio DNA sequencing rows a samplesheet date covers: every run up to the date plus
    READ_SET_WINDOW_DAYS. Returns (read set seq_date, [sequencing_id]). The read set's seq_date is
    the latest seq_date among them, or the samplesheet date when the OG has none (external reads),
    which v4.v_read_set_run_readiness then lists."""
    cur.execute(
        "SELECT s.sequencing_id, s.seq_date FROM v_sequencing_browse s "
        "WHERE s.og_id = %s AND s.technology IN ('PacBio HIFI', 'PacBio') AND s.library_type = 'pacbio' "
        "AND s.status IS DISTINCT FROM 'Failed' "  # a failed run has no reads (migration 069)
        "AND s.seq_date ~ '^[0-9]{6}$' "
        "AND to_date(s.seq_date, 'YYMMDD') <= to_date(%s, 'YYMMDD') + %s "
        "ORDER BY s.seq_date, s.sequencing_id",
        (og_id, sheet_date, READ_SET_WINDOW_DAYS))
    runs = cur.fetchall() or []
    return max((d for _, d in runs), default=sheet_date), [sid for sid, _ in runs]


def ensure_read_set(cur, og_id, sheet_date):
    """read_set_id for (og_id, samplesheet date), created with its runs if missing.

    An existing read set keeps the runs it has (genomes_db_review migration 093). They are the reads
    the assemblies on it were built from, and 093 corrected some from the k-mer evidence (OG838
    240813 is the 240813 cell alone, though the date rule also gives the 240520 one). Adding the
    date rule's runs on every push would undo that. A read set with no runs yet (one 093 or the
    ref_assembly trigger made before the OG's sequencing rows were loaded) gets them."""
    seq_date, runs = read_set_runs(cur, og_id, sheet_date)
    cur.execute(
        "SELECT rs.read_set_id, EXISTS (SELECT 1 FROM ref_read_set_run r WHERE r.read_set_id = rs.read_set_id) "
        "FROM ref_read_set rs WHERE rs.og_id = %s AND coalesce(rs.seq_date, '') = %s",
        (og_id, seq_date))
    found = cur.fetchone()
    if found and found[1]:
        return found[0]
    cur.execute(
        "INSERT INTO ref_read_set (og_id, seq_date) VALUES (%s, %s) "
        "ON CONFLICT (og_id, (coalesce(seq_date, ''))) DO UPDATE SET og_id = EXCLUDED.og_id "
        "RETURNING read_set_id",
        (og_id, seq_date))
    read_set_id = cur.fetchone()[0]
    for sid in runs:
        cur.execute("INSERT INTO ref_read_set_run (read_set_id, sequencing_id) VALUES (%s, %s) "
                    "ON CONFLICT DO NOTHING", (read_set_id, sid))
    return read_set_id


def upsert_read_set(cur, table, og_id, seq_date, values):
    """GenomeScope / HiFiAdapterFilt for one HiFi read set (v4.ref_read_set, migration 062).

    seq_date is the samplesheet's date (YYMMDD) and is required: a result whose reads are unknown is
    not written. All-None values remove the tool row, and the read set too once no tool row is left.
    Returns 'upserted', 'deleted' or 'unchanged'.
    """
    if table not in READ_SET_TABLES:
        raise ValueError(f"{table} is not a read-set table")
    unknown = set(values) - set(READ_SET_TABLES[table])
    if unknown:
        raise ValueError(f"{table} has no column(s) {sorted(unknown)}")
    og_id, seq_date = clean(og_id), clean(seq_date)
    if not og_id or not str(og_id).startswith("OG"):
        raise ValueError(f"bad og_id {og_id!r}")
    if seq_date is None:
        raise ValueError(f"no samplesheet date for {og_id}: pass the samplesheet so the read set is known")
    if not SEQ_DATE_RE.match(str(seq_date)):
        raise ValueError(f"bad seq_date {seq_date!r} for {og_id} (want YYMMDD)")
    seq_date = str(seq_date)
    vals = {c: clean(v) for c, v in values.items()}
    if all(v is None for v in vals.values()):
        key, _ = read_set_runs(cur, og_id, seq_date)
        cur.execute(f"DELETE FROM {table} t USING ref_read_set rs "
                    f"WHERE t.read_set_id = rs.read_set_id AND rs.og_id = %s AND rs.seq_date = %s",
                    (og_id, key))
        deleted = cur.rowcount
        # An assembly's read set stays (ref_assembly.read_set_id, migration 093), QC or not.
        cur.execute("DELETE FROM ref_read_set rs WHERE rs.og_id = %s AND rs.seq_date = %s "
                    "AND NOT EXISTS (SELECT 1 FROM ref_genomescope g WHERE g.read_set_id = rs.read_set_id) "
                    "AND NOT EXISTS (SELECT 1 FROM ref_hifi_adapter_filter f WHERE f.read_set_id = rs.read_set_id) "
                    "AND NOT EXISTS (SELECT 1 FROM ref_assembly a WHERE a.read_set_id = rs.read_set_id)",
                    (og_id, key))
        return "deleted" if deleted else "unchanged"
    read_set_id = ensure_read_set(cur, og_id, seq_date)
    cols = list(vals)
    cur.execute(
        f"INSERT INTO {table} (read_set_id, {', '.join(cols)}) "
        f"VALUES (%(read_set_id)s, {', '.join('%(' + c + ')s' for c in cols)}) "
        f"ON CONFLICT (read_set_id) DO UPDATE SET "
        f"{', '.join(f'{c} = EXCLUDED.{c}' for c in cols)}, loaded_at = now()",
        dict(vals, read_set_id=read_set_id))
    return "upserted"


# -- HiFiAdapterFilt per sequencing run (genomes_db_review migration 093) ---------------------------
#
# HiFiAdapterFilt filters each cell's BAM separately and writes one .stats file per BAM, so its
# count belongs to the sequencing row (ref_hifi_adapter_filter_run); the read-set figure in
# ref_hifi_adapter_filter is their sum. The rules below are the ones 093's backfill used
# (genomes_db_review/scripts/match_adapter_filter_runs.py).

STATS_MOVIE_RE = re.compile(r"(m\d+[a-z]?_(\d{6})_\d{6}(?:_s\d)?)")
STATS_TAG_RE = re.compile(r"\.(bc\d+(?:--bc\d+)?|unassigned)\.stats$")
STATS_COUNT_RE = re.compile(r"Number of adapter contaminated ccs reads:\s*(\d+)")


def parse_adapter_stats(text):
    """The adapter-contaminated read count in a HiFiAdapterFilt .stats file, or None."""
    m = STATS_COUNT_RE.search(text or "")
    return int(m.group(1)) if m else None


def adapter_stats_cell(path):
    """(movie, movie date YYMMDD, tag) from a .stats file name; tag is 'bcNNNN', 'unassigned' or
    'none' (unbarcoded). (None, None, tag) when the name carries no movie (old SSB_cell1_ccs.stats)."""
    name = Path(str(path)).name
    m = STATS_MOVIE_RE.search(name)
    t = STATS_TAG_RE.search(name)
    tag = t.group(1).split("--")[-1] if t else "none"
    return (m.group(1), m.group(2), tag) if m else (None, None, tag)


def pacbio_cells(cur, og_id):
    """The OG's non-failed PacBio DNA rows: [(sequencing_id, seq_date, barcode or None)]."""
    cur.execute(
        "SELECT DISTINCT s.sequencing_id, s.seq_date, q.barcode FROM v_sequencing_browse s "
        "LEFT JOIN hifi_read_qc q ON q.sequencing_id = s.sequencing_id "
        "WHERE s.og_id = %s AND s.technology IN ('PacBio HIFI', 'PacBio') AND s.library_type = 'pacbio' "
        "AND s.status IS DISTINCT FROM 'Failed' AND s.seq_date ~ '^[0-9]{6}$' "
        "ORDER BY s.seq_date, s.sequencing_id",
        (og_id,))
    return cur.fetchall() or []


def _days(a, b):
    from datetime import datetime
    return abs((datetime.strptime(a, "%y%m%d") - datetime.strptime(b, "%y%m%d")).days)


def _nearest(cells, movie_date):
    """The sequencing_id nearest the movie date within READ_SET_WINDOW_DAYS, or (None, why)."""
    best = sorted((_days(d, movie_date), sid) for sid, d, _ in cells if _days(d, movie_date) <= READ_SET_WINDOW_DAYS)
    if not best:
        return None, f"no sequencing row within {READ_SET_WINDOW_DAYS} days of the movie"
    ties = sorted({sid for n, sid in best if n == best[0][0]})
    if len(ties) > 1:
        return None, "tie: " + ", ".join(ties)
    return ties[0], ""


def map_adapter_stats(files, cells):
    """Assign one OG's .stats files to its sequencing rows.

    files: [(path, contam_reads)]. cells: pacbio_cells(). Returns ({sequencing_id: (contam_reads
    summed, [paths])}, [(path, why)] for files left out).
    - Copies of one file (hifi/ and hifi/qc_stats/, or several assembly folders) count once, and
      copies that disagree raise ValueError.
    - A barcoded file goes to the row with that barcode (hifi_read_qc) nearest the movie date within
      READ_SET_WINDOW_DAYS; when no row has the barcode, to the nearest row by date.
    - An unassigned or unbarcoded file is part of the same cell: it follows the movie's barcoded
      file, else the nearest row by date. One whose movie's barcoded files went to several rows is
      left out, silently when its count is 0 (it adds nothing to either).
    """
    seen = {}
    left_out = []
    for path, n in files:
        movie, movie_date, tag = adapter_stats_cell(path)
        if movie is None:
            left_out.append((path, "no PacBio movie in the file name"))
            continue
        if n is None:
            left_out.append((path, "no adapter count in the file"))
            continue
        key = (movie, tag)
        if key in seen and seen[key][1] != n:
            raise ValueError(f"copies of {movie} {tag} disagree: {seen[key][0]} has {seen[key][1]}, {path} has {n}")
        seen.setdefault(key, (path, n, movie_date))

    by_movie = {}
    for (movie, tag), f in seen.items():
        by_movie.setdefault(movie, {})[tag] = f

    assigned = []
    for movie, tags in sorted(by_movie.items()):
        barcoded = set()
        for tag, (path, n, movie_date) in sorted(tags.items()):
            if not tag.startswith("bc"):
                continue
            with_bc = [c for c in cells if c[2] == tag]
            sid, why = _nearest(with_bc or cells, movie_date)
            if sid:
                barcoded.add(sid)
                assigned.append((sid, path, n))
            else:
                left_out.append((path, why))
        for tag, (path, n, movie_date) in sorted(tags.items()):
            if tag.startswith("bc"):
                continue
            if len(barcoded) == 1:
                assigned.append((next(iter(barcoded)), path, n))
            elif len(barcoded) > 1:
                if n:
                    left_out.append((path, "the movie's barcoded files went to several rows"))
            else:
                sid, why = _nearest(cells, movie_date)
                if sid:
                    assigned.append((sid, path, n))
                else:
                    left_out.append((path, why))

    per_run = {}
    for sid, path, n in assigned:
        total, paths = per_run.get(sid, (0, []))
        per_run[sid] = (total + n, paths + [path])
    return per_run, left_out


def upsert_adapter_filter_run(cur, sequencing_id, contam_reads, source=None):
    """One sequencing row's HiFiAdapterFilt count (ref_hifi_adapter_filter_run). A cell's count does
    not change when more data arrives, so a re-push replaces it with the same figure."""
    n = parse_int(contam_reads)
    if n is None or n < 0:
        raise ValueError(f"bad contam_reads {contam_reads!r} for {sequencing_id}")
    cur.execute(
        "INSERT INTO ref_hifi_adapter_filter_run (sequencing_id, contam_reads, source) VALUES (%s, %s, %s) "
        "ON CONFLICT (sequencing_id) DO UPDATE SET contam_reads = EXCLUDED.contam_reads, "
        "source = EXCLUDED.source, loaded_at = now()",
        (sequencing_id, n, source))
    return "upserted"


def write_adapter_filter(db_params, read_set_rows, stats_files):
    """HiFiAdapterFilt in one transaction: per run, then per read set.

    read_set_rows: as write_read_set_rows() for ref_hifi_adapter_filter.
    stats_files: {og_id: [(path, contam_reads)]}, every .stats file found for the OG.
    A file that maps to no sequencing row is reported and left out, not fatal:
    v_ref_read_set_readiness lists the runs still without a count.
    """
    rows = list(read_set_rows)

    def writes(cur):
        counts = {}
        for og_id, files in sorted(stats_files.items()):
            per_run, left_out = map_adapter_stats(files, pacbio_cells(cur, og_id))
            for sid, (n, paths) in sorted(per_run.items()):
                upsert_adapter_filter_run(cur, sid, n, "; ".join(sorted(paths)))
                counts["run upserted"] = counts.get("run upserted", 0) + 1
            for path, why in left_out:
                print(f"⚠️  {og_id}: {path} left out: {why}")
                counts["file left out"] = counts.get("file left out", 0) + 1
        for og_id, seq_date, values in rows:
            outcome = upsert_read_set(cur, "ref_hifi_adapter_filter", og_id, seq_date, values)
            counts[outcome] = counts.get(outcome, 0) + 1
        return counts
    return write_all(db_params, writes, "ref_hifi_adapter_filter_run + ref_hifi_adapter_filter")


# -- Hi-C runs of an assembly (genomes_db_review migration 094) -------------------------------------
# An Illumina read name starts instrument:run_number:flowcell:lane, and bwa mem and pairtools keep the
# names, so the flowcells among an assembly's Hi-C BAM records are the runs its Hi-C came from. The
# database maps a flowcell to its run (sequencing_run.flowcell) and derives the links
# (set_ref_assembly_hic_flowcells, v_ref_assembly_hic_rule). Only the standard library is used here,
# so this runs in the psycopg2 container, which has no samtools.

BAM_SAMPLE_BYTES = 4_000_000  # ~50,000 records of a coordinate-sorted BAM: every input run shows


def bam_flowcells(path, max_bytes=BAM_SAMPLE_BYTES):
    """{flowcell: records} over the first max_bytes of a BAM. Raises ValueError if it is not a BAM
    or no read name carries a flowcell (a BAM whose names were rewritten cannot say)."""
    import struct
    import zlib
    with open(Path(path).expanduser(), "rb") as fh:
        data = fh.read(max_bytes)
    out = bytearray()
    i = 0
    while i + 18 <= len(data) and data[i:i + 4] == b"\x1f\x8b\x08\x04":
        bsize = struct.unpack("<H", data[i + 16:i + 18])[0] + 1
        if i + bsize > len(data):
            break
        out += zlib.decompress(data[i + 18:i + bsize - 8], -15)
        i += bsize
    if out[:4] != b"BAM\x01":
        raise ValueError(f"{path}: not a BAM file")
    p = 8 + struct.unpack("<i", out[4:8])[0]
    n_ref = struct.unpack("<i", out[p:p + 4])[0]
    p += 4
    for _ in range(n_ref):
        p += 4 + struct.unpack("<i", out[p:p + 4])[0] + 4
    counts = {}
    while p + 4 <= len(out):
        bs = struct.unpack("<i", out[p:p + 4])[0]
        if p + 4 + bs > len(out):
            break
        l_name = out[p + 12]
        f = bytes(out[p + 36:p + 36 + l_name - 1]).decode(errors="replace").split(":")
        if len(f) >= 7:
            counts[f[2]] = counts.get(f[2], 0) + 1
        p += 4 + bs
    if not counts:
        raise ValueError(f"{path}: no Illumina read names with a flowcell in the first {max_bytes:,} bytes")
    return counts


def set_hic_flowcells(cur, og_id, seq_date, version, flowcells):
    """Record an assembly's Hi-C flowcells ({flowcell: records}); the database derives its Hi-C links.
    Returns the problems the database reports (a flowcell no run has; a run with no non-failed Hi-C
    row of the OG): until they are fixed the assembly has no Hi-C link."""
    cur.execute("SELECT assembly_id FROM ref_assembly WHERE og_id = %s AND seq_date = %s AND version = %s",
                (og_id, str(seq_date), version))
    row = cur.fetchone()
    if not row:
        raise ValueError(f"no ref_assembly {og_id} {seq_date}.{version}: push the assembly's results first")
    fcs = sorted(flowcells)
    cur.execute("SELECT set_ref_assembly_hic_flowcells(%s, %s, %s)", (row[0], fcs, [flowcells[f] for f in fcs]))
    return list((cur.fetchone() or [[]])[0] or [])


def write_hic_flowcells(db_params, items):
    """items: [((og_id, seq_date, version), {flowcell: records})]. One transaction, committed, then
    raises RuntimeError listing every problem, so the caller exits non-zero while the resolvable
    assemblies stay linked (v_ref_assembly_hic_readiness lists the rest)."""
    items = list(items)
    problems = []

    def writes(cur):
        counts = {}
        for (og_id, seq_date, version), flowcells in items:
            found = set_hic_flowcells(cur, og_id, seq_date, version, flowcells)
            problems.extend(f"{og_id} {seq_date}.{version}: {p}" for p in found)
            outcome = "unresolved" if found else "linked"
            counts[outcome] = counts.get(outcome, 0) + 1
        return counts
    write_all(db_params, writes, "ref_assembly_hic_flowcell")
    if problems:
        raise RuntimeError("Hi-C runs not linked:\n  " + "\n  ".join(problems))


def upsert_hifi_read_qc(cur, row):
    """One PacBio HiFi library-on-run row. Links the one sequencing row of that library on that run,
    rebuilding the tube as og_id||tissue||'_'||ext_type||'_'||lib_code (NULL if none or several)."""
    r = {c: clean(row.get(c)) for c in HIFI_READ_QC_KEY + HIFI_READ_QC_COLS}
    blank = [c for c in HIFI_READ_QC_KEY if not r[c]]
    if blank:
        raise ValueError(f"hifi_read_qc row with blank {blank}: {r}")
    cur.execute(
        f"INSERT INTO hifi_read_qc ({', '.join(HIFI_READ_QC_KEY)}, sequencing_id, {', '.join(HIFI_READ_QC_COLS)}) "
        f"SELECT {', '.join('%(' + c + ')s' for c in HIFI_READ_QC_KEY)}, "
        "(SELECT CASE WHEN count(*) = 1 THEN min(s.sequencing_id) END FROM sequencing s "
        " WHERE s.library_tube_id = %(og_id)s || %(tissue)s || '_' || %(ext_type)s || '_' || %(lib_code)s "
        "   AND s.run_id = %(run_id)s), "
        f"{', '.join('%(' + c + ')s' for c in HIFI_READ_QC_COLS)} "
        f"ON CONFLICT ({', '.join(HIFI_READ_QC_KEY)}) DO UPDATE SET sequencing_id = EXCLUDED.sequencing_id, "
        f"{', '.join(f'{c} = EXCLUDED.{c}' for c in HIFI_READ_QC_COLS)}, loaded_at = now()",
        r)
    return "upserted"


# -- transactions -----------------------------------------------------------------------------------

def connect(db_params):
    import psycopg2
    return psycopg2.connect(**db_params)


def write_all(db_params, writes, label):
    """Run writes(cur) -> {outcome: count} in one transaction. Commit only if it all succeeds."""
    conn = connect(db_params)
    try:
        with conn:
            with conn.cursor() as cur:
                counts = writes(cur)
    finally:
        conn.close()
    summary = ", ".join(f"{n} {k}" for k, n in sorted(counts.items())) or "nothing to write"
    print(f"✅ {label}: {summary} (one transaction, committed)")
    return counts


def write_tool_rows(db_params, table, rows):
    """rows: iterable of (key dict, values dict). One transaction for the whole file."""
    rows = list(rows)

    def writes(cur):
        counts = {}
        for key, values in rows:
            outcome = upsert_metrics(cur, table, key, values)
            counts[outcome] = counts.get(outcome, 0) + 1
        return counts
    return write_all(db_params, writes, table)


def write_read_set_rows(db_params, table, rows):
    """rows: iterable of (og_id, samplesheet date YYMMDD, values dict). One transaction; a row with
    no date fails the whole file."""
    rows = list(rows)

    def writes(cur):
        counts = {}
        for og_id, seq_date, values in rows:
            outcome = upsert_read_set(cur, table, og_id, seq_date, values)
            counts[outcome] = counts.get(outcome, 0) + 1
        return counts
    return write_all(db_params, writes, table)


def write_hifi_read_qc_rows(db_params, rows):
    rows = list(rows)

    def writes(cur):
        for row in rows:
            upsert_hifi_read_qc(cur, row)
        return {"upserted": len(rows)}
    return write_all(db_params, writes, "hifi_read_qc")


def main_guard(main):
    """Run main(); on any error print it and exit 1, so pipelines see the failure."""
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001 - report everything, then fail
        print(f"❌ Error: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(1)
