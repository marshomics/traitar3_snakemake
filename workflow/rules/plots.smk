# QC / summary / biological figures. All rules use the SEPARATE plots.yaml env
# (not traitar3.yaml), so adding them does not change the annotation/prediction
# software-env signatures and cannot re-trigger those completed jobs.

# --- annotation QC statistics (batched, from existing count files) -----------

rule qc_stats_batch:
    """Per-genome annotation stats for one batch (reads the small count files)."""
    input:
        counts=lambda wc: [counts_path(s) for s in BATCHES[int(wc.batch)]],
        manifest=f"{OUT}/predict/_manifests/batch_{{batch}}.counts.txt",
        accs=f"{OUT}/db/model_pfam_accessions.txt",
    output:
        stats=f"{OUT}/qc/batch_{{batch}}.stats.tsv",
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb=res("qc_stats", "mem_mb"),
        runtime=res("qc_stats", "runtime"),
    shell:
        "python workflow/scripts/qc_stats.py --counts-list {input.manifest} "
        "--model-accessions {input.accs} --out {output.stats}"


rule qc_stats_gather:
    """Concatenate per-batch QC stats into one genome x metric table."""
    input:
        parts=[f"{OUT}/qc/batch_{b}.stats.tsv" for b in range(N_BATCHES)],
    output:
        table=f"{OUT}/qc/annotation_stats.tsv",
    run:
        os.makedirs(os.path.dirname(output.table), exist_ok=True)
        with open(output.table, "w") as out:
            for i, p in enumerate(input.parts):
                with open(p) as fh:
                    header = fh.readline()
                    if i == 0:
                        out.write(header)
                    for line in fh:
                        out.write(line)


# --- figures -----------------------------------------------------------------

rule plot_qc:
    input:
        stats=f"{OUT}/qc/annotation_stats.tsv",
        combined=f"{OUT}/predictions_majority-vote_combined.tsv",
    output:
        figs=expand(f"{PLOTS_DIR}/{{fig}}.{{fmt}}", fig=QC_FIGS, fmt=PLOT_FORMATS),
        caption=f"{PLOTS_DIR}/FIGURES.qc.md",
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb=res("plots", "mem_mb"),
        runtime=res("plots", "runtime"),
    params:
        fmts=",".join(PLOT_FORMATS),
        dpi=PLOT_DPI,
        low=PLOT_LOW_FAM,
        outdir=PLOTS_DIR,
    shell:
        "python workflow/scripts/plot_qc.py --annotation-stats {input.stats} "
        "--combined {input.combined} --outdir {params.outdir} "
        "--formats {params.fmts} --dpi {params.dpi} --low-families {params.low}"


rule plot_summary:
    input:
        combined=f"{OUT}/predictions_majority-vote_combined.tsv",
        pp_single=f"{OUT}/predict/{PRIMARY['token']}.single_votes.tsv",
        pgl_single=f"{OUT}/predict/{SECONDARY['token']}.single_votes.tsv",
        model=PRIMARY_ARCHIVE,
    output:
        figs=expand(f"{PLOTS_DIR}/{{fig}}.{{fmt}}", fig=SUMMARY_FIGS, fmt=PLOT_FORMATS),
        caption=f"{PLOTS_DIR}/FIGURES.summary.md",
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb=res("plots", "mem_mb"),
        runtime=res("plots", "runtime"),
    params:
        fmts=",".join(PLOT_FORMATS),
        dpi=PLOT_DPI,
        outdir=PLOTS_DIR,
        thr=config["predict"]["majority_threshold"],
    shell:
        "python workflow/scripts/plot_summary.py --combined {input.combined} "
        "--phypat-single {input.pp_single} --phypatpgl-single {input.pgl_single} "
        "--model-tar {input.model} --outdir {params.outdir} "
        "--formats {params.fmts} --dpi {params.dpi} --majority-threshold {params.thr}"


rule plot_biological:
    input:
        combined=f"{OUT}/predictions_majority-vote_combined.tsv",
        model=PRIMARY_ARCHIVE,
    output:
        figs=expand(f"{PLOTS_DIR}/{{fig}}.{{fmt}}", fig=BIO_FIGS, fmt=PLOT_FORMATS),
        caption=f"{PLOTS_DIR}/FIGURES.biological.md",
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb=res("plots", "mem_mb"),
        runtime=res("plots", "runtime"),
    params:
        fmts=",".join(PLOT_FORMATS),
        dpi=PLOT_DPI,
        outdir=PLOTS_DIR,
    shell:
        "python workflow/scripts/plot_biological.py --combined {input.combined} "
        "--model-tar {input.model} --outdir {params.outdir} "
        "--formats {params.fmts} --dpi {params.dpi}"
