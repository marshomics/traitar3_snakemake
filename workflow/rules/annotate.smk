# Per-genome annotation: hmmsearch against (subset) Pfam-A, then collapse to a
# per-genome Pfam count vector. One job per genome = maximal parallelism.
#
# The raw domtblout is a throwaway, so it's written to the node-local job temp
# (resources.tmpdir, = $TMPDIR on the cluster) rather than into the shared output
# directory. That avoids depending on the shared filesystem having created/
# propagated results/counts/ before the job starts (an NFS race that otherwise
# makes a fraction of jobs fail with "Failed to open ... for writing"), and it's
# faster. The job also creates its own output directory to be safe.

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
    shell:
        r"""
        mkdir -p "$(dirname {output.counts})"
        domtbl="{resources.tmpdir}/{wildcards.sample}.domtbl.tmp"
        hmmsearch {params.extra} --cpu {threads} \
            --domtblout "$domtbl" {input.hmm} {input.faa} > /dev/null
        python workflow/scripts/hmmer_to_counts.py "$domtbl" {output.counts} \
            --evalue-max {params.evalue_max} --score-min {params.score_min}
        rm -f "$domtbl"
        """
