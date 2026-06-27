# Per-genome annotation: hmmsearch against (subset) Pfam-A, then collapse to a
# per-genome Pfam count vector. One job per genome = maximal parallelism. The
# raw domtblout is written to a job-local temp file and deleted in the same job
# so 300k genomes don't leave 300k large intermediates behind.

rule annotate:
    input:
        faa=lambda wc: SAMPLE2FAA[wc.sample],
        hmm=HMM_DB,
    output:
        counts=f"{OUT}/counts/{{sample}}.pfam_counts.tsv",
    conda:
        "../envs/traitar3.yaml"
    threads: config["hmmsearch"]["threads"]
    resources:
        mem_mb=res("annotate", "mem_mb"),
        runtime=res("annotate", "runtime"),
    params:
        extra=config["hmmsearch"]["extra"],
        evalue_max=config["filter"]["evalue_max"],
        score_min=config["filter"]["score_min"],
        domtbl=lambda wc, output: output.counts + ".domtbl.tmp",
    shell:
        r"""
        hmmsearch {params.extra} --cpu {threads} \
            --domtblout {params.domtbl} {input.hmm} {input.faa} > /dev/null
        python workflow/scripts/hmmer_to_counts.py {params.domtbl} {output.counts} \
            --evalue-max {params.evalue_max} --score-min {params.score_min}
        rm -f {params.domtbl}
        """
