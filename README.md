# traitar3-smk

A Snakemake reimplementation of the [traitar3](https://github.com/nick-youngblut/traitar3)
`from_genes` phenotype pipeline, built to run hundreds of thousands of genomes
on a cluster. Input is amino-acid FASTA (`*.faa`) already called by Prodigal (or
any gene caller); output is traitar's 67-trait phenotype calls.

Per genome the steps are: `hmmsearch` against Pfam-A, filter to the best domain
hit per (protein, Pfam), reduce to a per-genome Pfam count vector, then score the
phypat and phypat+PGL committees and merge them into combined calls. Each genome
is independent, so annotation is one job per genome; prediction runs on batches
so a model is loaded once per batch instead of once per genome. This replaces
traitar3's single-process design, whose in-memory `genomes x ~14,800-Pfam`
matrix and per-cell pandas aggregation don't survive 300k genomes.

## Why this is faithful (and where it deviates on purpose)

The annotation, counting, prediction math, and the phypat/phypat+PGL merge are
ported directly from traitar3's own code. Validation against the repo's bundled
reference data (`test/run_offline_check.sh`):

- per-genome Pfam counts reproduce traitar's `summary.dat` exactly;
- per-model votes reproduce the decision values in traitar's committed
  `predictions_raw.txt` exactly;
- the combined majority-vote calls reproduce `predictions_majority-vote_combined.txt`
  exactly.

One deliberate deviation: traitar3's literal Python-3 code sets the majority cutoff
to `votes >= k/2 + 1`, which for `k=5` is `>= 3.5` (4 of 5 voters) — a regression
introduced when integer division (`/`) changed meaning from Python 2. The original
traitar and traitar3's own reference outputs use a true majority of `>= 3`. The
default here is `>= 3` (`config: predict.majority_threshold`); set it to `3.5`
(or pass `--literal-traitar3` to the predict script) to reproduce traitar3's
literal numbers.

## Layout

```
config/config.yaml          # all settings
workflow/Snakefile          # entry point
workflow/rules/             # common / download / annotate / predict
workflow/scripts/           # the ported logic (counts, predict, merge, subset)
workflow/envs/traitar3.yaml # conda env used by every rule (--use-conda)
profiles/sge/               # SGE/UGE cluster profile (qsub submit + status)
resources/models/           # phypat.tar.gz, phypat+PGL.tar.gz (bundled)
test/                       # 2 sample genomes + reference outputs + check
```

## Install

The controller environment needs Snakemake (>=8) and, for the cluster, the
generic executor plugin:

```bash
conda create -n smk -c conda-forge -c bioconda \
    "snakemake>=8" snakemake-executor-plugin-cluster-generic
conda activate smk
```

Each rule pulls its own tools (HMMER, pandas, numpy) from
`workflow/envs/traitar3.yaml` when you pass `--use-conda`, so nothing else needs
installing by hand.

## Pre-build the conda environment (before submitting to the cluster)

Compute nodes often have no internet, so build the environments once on the login
node and let the jobs reuse them. The workflow defines two: `traitar3.yaml`
(HMMER, numpy, pandas) for annotation and prediction, and `plots.yaml`
(matplotlib, seaborn, scipy) for the figures. `--conda-create-envs-only` builds
both.

```bash
# on the login node, from the repo root
snakemake --use-conda --conda-create-envs-only --conda-frontend mamba --cores 4
```

This creates the env under `.snakemake/conda/` (a content hash of the YAML). Jobs
find it by that same hash, so as long as you submit from the same repo directory
on the shared filesystem, the compute nodes reuse it and never touch the network.
To keep envs in a stable shared location instead, add the same
`--conda-prefix /shared/path/conda-envs` to both this command and the run.

While you're on the login node, also pre-stage the database (the `download_pfam`
rule needs internet too):

```bash
snakemake --use-conda --cores 4 results/db/Pfam-A.model-subset.hmm
```

Then submit; everything is already built, so jobs run offline:

```bash
snakemake --use-conda --workflow-profile profiles/sge
```

## Check it works (no database needed)

```bash
bash test/run_offline_check.sh        # prints PASSED
```

## Configure

Edit `config/config.yaml`:

- `input_dir` / `faa_extension` — directory scanned for `*.faa`. The sample name
  is the filename without the extension. For genomes spread across directories,
  set `samples_tsv` to a `sample<TAB>faa` table instead.
- `pfam.subset: true` — search only the ~3,400 Pfam families that carry non-zero
  predictor weight (of 14,831 in Pfam 27.0). Predictions are identical, since a
  zero-weight family can't change a linear decision, and `hmmsearch` scans ~4x
  fewer profiles. Set `false` to search the full database.
- `batch_size` — genomes per prediction job (1000–5000 is sensible).
- `resources:` — per-rule memory (MB) and runtime (minutes) for the cluster.

## Run

```bash
# local
snakemake --use-conda --cores 16

# SGE / UGE cluster (Pfam download + subset happen automatically on first run).
# The command differs by Snakemake version (check `snakemake --version`):
snakemake --use-conda --profile profiles/sge            # Snakemake 7
snakemake --use-conda --workflow-profile profiles/sge   # Snakemake 8+
#   override the parallel-environment name if yours isn't "smp":
SGE_PE=threads snakemake --use-conda --profile profiles/sge
#   if mamba isn't on the cluster, add: --conda-frontend conda
```

The `profiles/sge/` directory holds both: `config.yaml` (Snakemake 7, built-in
`--cluster`) and `config.v8+.yaml` (Snakemake 8, the cluster-generic executor
plugin, which needs `pip install snakemake-executor-plugin-cluster-generic`).
Snakemake 8 picks the `v8+` file automatically. On v7 use `--profile`, not
`--workflow-profile` (v7 mishandles profile script paths with the latter).

First run downloads Pfam-A 27.0 (~1.2 GB) into `resources/pfam/` and builds the
model subset once; both are reused thereafter.

## Outputs (`results/`)

- `predictions_majority-vote_combined.tsv` — genomes x 67 traits, encoded
  `0` negative, `1` phypat only, `2` phypat+PGL only, `3` both.
- `predictions_single-votes_combined.tsv` — summed vote counts (0–10).
- `predictions_flat_*_combined.tsv` — long-format, one non-zero call per row.
- `predictions_reliable.tsv` — genomes x 67, `1` where phypat+PGL called the trait
  positive (combined code 2 or 3), the paper's "reliable" set; `predictions_reliable_flat.tsv`
  is the tidy list of those calls. Add more sets (e.g. `strict: [3]`) under
  `derived_call_sets` in `config.yaml`.
- `predict/phypat.{single_votes,majority_vote}.tsv` and the `phypatPGL`
  equivalents — per-model tables.
- `counts/<sample>.pfam_counts.tsv` — per-genome Pfam presence (sparse).
- `plots/*.png` and `plots/*.svg` — publication figures (see below); SVG text is
  editable in Illustrator/Inkscape. `qc/annotation_stats.tsv` holds the
  per-genome QC metrics behind them.

## Figures

`results/plots/` holds nine figures, each as a PNG (600 DPI) + editable SVG pair,
with auto-written captions in `plots/FIGURES.*.md`. Every figure aggregates
across genomes (histograms, ECDFs, 2D density, per-trait bars, a 67×67
trait-correlation heatmap), so they stay readable and cheap at 340k proteomes —
each of the three plot jobs runs in seconds using under 1 GB of RAM.

QC (annotation quality): `qc_annotation_completeness` (Pfam families, depth and
model-feature coverage per genome), `qc_completeness_vs_calls` (density of
coverage vs number of calls, to check calls aren't driven by sparse genomes),
`qc_calls_per_genome` (positive traits per genome, with a zero-call flag).

Summary (prediction overview): `summary_trait_prevalence` (per-trait prevalence,
phypat vs phypat+PGL), `summary_model_agreement` (per-trait 0/1/2/3 concordance),
`summary_confidence` (committee vote support overall and among positive calls).

Biological: `bio_prevalence_by_category` (prevalence grouped by trait category),
`bio_trait_cooccurrence` (clustered 67×67 co-occurrence heatmap), and
`bio_key_traits_composition` (oxygen relationship, Gram stain, morphology).

Toggle or tune via the `plots:` block in `config/config.yaml` (`enabled`,
`formats`, `dpi`). The plot rules use a separate conda env
(`workflow/envs/plots.yaml`) so they never alter the annotation env or re-trigger
those jobs.

## Scaling to ~300k genomes

The design is one annotation job per genome (embarrassingly parallel) plus
`ceil(N / batch_size)` prediction jobs. Practical notes:

- **Throughput.** `hmmsearch` dominates. With `pfam.subset: true` a ~3 Mb
  proteome is a few minutes; against full Pfam-A budget ~10–15 min. 300k genomes
  is ~25k–50k CPU-hours — about a day on ~1,000 cores.
- **Scheduler load.** Cap concurrency with `jobs:` in the profile and keep
  `max-jobs-per-second` modest. Building a 300k-job DAG takes a few minutes; let
  it run. Consider splitting into a few invocations by `input_dir` subfolders if
  your scheduler dislikes one giant submission.
- **Storage / inodes.** Raw `domtblout` is written to a job-local temp file and
  deleted in the same job, so only the small `counts/` vectors persist (one tiny
  TSV per genome). The final tables are `genomes x 67`.
- **Restartability.** Re-running resumes from completed outputs;
  `restart-times: 2` retries transient cluster failures and `keep-going: true`
  lets the run finish around a single bad genome.
