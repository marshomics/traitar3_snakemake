# Sample discovery, model bookkeeping and shared helpers.
import os
import re
import glob
import csv
import math


def _sanitize(name: str) -> str:
    """traitar-style sample-name sanitisation: keep [A-Za-z0-9._-], replace rest."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def _discover_samples(cfg):
    """Return {sample_name: faa_path}, from samples_tsv if given, else by glob."""
    mapping = {}
    tsv = cfg.get("samples_tsv") or ""
    if tsv:
        with open(tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None or "sample" not in reader.fieldnames or "faa" not in reader.fieldnames:
                raise WorkflowError("samples_tsv must have a header with columns: sample<TAB>faa")
            for row in reader:
                s = _sanitize(row["sample"].strip())
                _add_sample(mapping, s, row["faa"].strip())
    else:
        ext = cfg["faa_extension"]
        pattern = os.path.join(cfg["input_dir"], "*" + ext)
        for path in sorted(glob.glob(pattern)):
            s = _sanitize(os.path.basename(path)[: -len(ext)] if path.endswith(ext)
                          else os.path.splitext(os.path.basename(path))[0])
            _add_sample(mapping, s, path)
    if not mapping:
        raise WorkflowError(
            "No input genomes found. Set config['input_dir']/faa_extension (or "
            "samples_tsv) so that at least one .faa file is discovered."
        )
    return mapping


def _add_sample(mapping, sample, faa):
    if sample in mapping and mapping[sample] != faa:
        raise WorkflowError(
            f"Two input files map to the same sample name '{sample}': "
            f"{mapping[sample]} and {faa}. Rename one or use samples_tsv."
        )
    mapping[sample] = faa


# --- resolve config ----------------------------------------------------------
OUT = config["outdir"]

SAMPLE2FAA = _discover_samples(config)
SAMPLES = sorted(SAMPLE2FAA)
N_SAMPLES = len(SAMPLES)

BATCH_SIZE = int(config["batch_size"])
N_BATCHES = max(1, math.ceil(N_SAMPLES / BATCH_SIZE))
BATCHES = {b: SAMPLES[b * BATCH_SIZE:(b + 1) * BATCH_SIZE] for b in range(N_BATCHES)}

MODELS = config["models"]
TOKENS = [m["token"] for m in MODELS]
TOKEN2ARCHIVE = {m["token"]: m["archive"] for m in MODELS}
TOKEN2NAME = {m["token"]: m["name"] for m in MODELS}
PRIMARY = next((m for m in MODELS if m.get("role") == "primary"), MODELS[0])
SECONDARY = next((m for m in MODELS if m.get("role") == "secondary"), None)
HAS_SECONDARY = SECONDARY is not None

# database the annotation step searches against (full Pfam-A or model subset)
HMM_DB = config["pfam"]["subset_hmm"] if config["pfam"]["subset"] else config["pfam"]["hmm"]

VOTE_TYPES = ["single_votes", "majority_vote"]

# write full long-format flat files, or header-only placeholders (they are huge)
WRITE_FLAT = bool(config.get("write_flat", True))

# --- figures -----------------------------------------------------------------
_PLOTS_CFG = config.get("plots", {})
PLOTS_ENABLED = bool(_PLOTS_CFG.get("enabled", True)) and HAS_SECONDARY
PLOT_FORMATS = list(_PLOTS_CFG.get("formats", ["png", "svg"]))
PLOT_DPI = int(_PLOTS_CFG.get("dpi", 600))
PLOT_LOW_FAM = int(_PLOTS_CFG.get("low_families", 50))
PLOTS_DIR = f"{OUT}/plots"
QC_FIGS = ["qc_annotation_completeness", "qc_completeness_vs_calls", "qc_calls_per_genome"]
SUMMARY_FIGS = ["summary_trait_prevalence", "summary_model_agreement", "summary_confidence"]
BIO_FIGS = ["bio_prevalence_by_category", "bio_trait_cooccurrence", "bio_key_traits_composition"]
ALL_FIGS = QC_FIGS + SUMMARY_FIGS + BIO_FIGS
PRIMARY_ARCHIVE = TOKEN2ARCHIVE[PRIMARY["token"]]

# named binary call sets derived from the combined 0/1/2/3 matrix (needs both models)
DERIVED = dict(config.get("derived_call_sets", {})) if HAS_SECONDARY else {}


def counts_path(sample):
    return f"{OUT}/counts/{sample}.pfam_counts.tsv"


def res(rule_name, key):
    """Look up a per-rule resource value from config['resources']."""
    return config["resources"][rule_name][key]


def final_targets():
    targets = []
    # per-model gathered tables (always produced)
    for tok in TOKENS:
        for vt in VOTE_TYPES:
            targets.append(f"{OUT}/predict/{tok}.{vt}.tsv")
    # combined calls when a secondary model is configured
    if HAS_SECONDARY:
        targets += [
            f"{OUT}/predictions_majority-vote_combined.tsv",
            f"{OUT}/predictions_single-votes_combined.tsv",
            f"{OUT}/predictions_flat_majority-votes_combined.tsv",
            f"{OUT}/predictions_flat_single-votes_combined.tsv",
        ]
        # derived binary call sets (e.g. the phypat+PGL "reliable" set)
        for name in DERIVED:
            targets.append(f"{OUT}/predictions_{name}.tsv")
            targets.append(f"{OUT}/predictions_{name}_flat.tsv")
    # QC / summary / biological figures
    if PLOTS_ENABLED:
        for fig in ALL_FIGS:
            for fmt in PLOT_FORMATS:
                targets.append(f"{PLOTS_DIR}/{fig}.{fmt}")
    return targets


wildcard_constraints:
    sample=r"[A-Za-z0-9._-]+",
    batch=r"\d+",
    token=r"[A-Za-z0-9]+",
    vt=r"single_votes|majority_vote",
    name=r"[A-Za-z0-9_]+",
