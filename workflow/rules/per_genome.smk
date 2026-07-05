# Optional per-genome majority-vote output: one file per input genome, written
# from the per-batch prediction tables that already exist. Purely additive -- it
# reads predict_batch outputs and does not change merge, plots, or any other
# stage. A per-batch sentinel is the tracked output so the workflow doesn't have
# to declare 350k individual files.

rule per_genome_batch:
    input:
        primary=f"{OUT}/predict/{PRIMARY['token']}/batch_{{batch}}.majority_vote.tsv",
        secondary=f"{OUT}/predict/{SECONDARY['token']}/batch_{{batch}}.majority_vote.tsv",
        model=PRIMARY_ARCHIVE,
    output:
        marker=f"{OUT}/per_genome_markers/batch_{{batch}}.done",
    conda:
        "../envs/traitar3.yaml"
    resources:
        mem_mb=res("per_genome", "mem_mb"),
        runtime=res("per_genome", "runtime"),
    params:
        outdir=lambda wc: (f"{OUT}/per_genome/batch_{wc.batch}"
                           if PER_GENOME_SHARD else f"{OUT}/per_genome"),
        pname=PRIMARY["name"],
        sname=SECONDARY["name"],
    shell:
        "python workflow/scripts/split_per_genome.py "
        "--primary {input.primary} --secondary {input.secondary} "
        "--model-tar {input.model} "
        "--primary-name '{params.pname}' --secondary-name '{params.sname}' "
        "--outdir {params.outdir} --marker {output.marker}"
