# Batched prediction, gather, and phypat/phypat+PGL merge.

# --- manifests (tiny, run on the submit host) --------------------------------

rule batch_counts_manifest:
    """List the per-genome count files belonging to one batch."""
    output:
        manifest=f"{OUT}/predict/_manifests/batch_{{batch}}.counts.txt",
    run:
        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
        with open(output.manifest, "w") as fh:
            for s in BATCHES[int(wildcards.batch)]:
                fh.write(counts_path(s) + "\n")


rule gather_manifest:
    """List the per-batch prediction tables for one model + vote type."""
    output:
        manifest=f"{OUT}/predict/_manifests/{{token}}.{{vt}}.gather.txt",
    run:
        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
        with open(output.manifest, "w") as fh:
            for b in range(N_BATCHES):
                fh.write(f"{OUT}/predict/{wildcards.token}/batch_{b}.{wildcards.vt}.tsv\n")


# --- prediction --------------------------------------------------------------

rule predict_batch:
    """Score one batch of genomes with one phenotype model."""
    input:
        counts=lambda wc: [counts_path(s) for s in BATCHES[int(wc.batch)]],
        manifest=f"{OUT}/predict/_manifests/batch_{{batch}}.counts.txt",
        archive=lambda wc: TOKEN2ARCHIVE[wc.token],
    output:
        sv=f"{OUT}/predict/{{token}}/batch_{{batch}}.single_votes.tsv",
        mv=f"{OUT}/predict/{{token}}/batch_{{batch}}.majority_vote.tsv",
    conda:
        "../envs/traitar3.yaml"
    resources:
        mem_mb=res("predict", "mem_mb"),
        runtime=res("predict", "runtime"),
    params:
        prefix=lambda wc, output: output.sv[: -len(".single_votes.tsv")],
        k=config["predict"]["voters"],
        thr=config["predict"]["majority_threshold"],
    shell:
        "python workflow/scripts/predict_models.py --model {input.archive} "
        "--counts-list {input.manifest} --out-prefix {params.prefix} "
        "--voters {params.k} --majority-threshold {params.thr}"


rule gather:
    """Concatenate per-batch tables into one genome x phenotype table per model."""
    input:
        parts=lambda wc: [f"{OUT}/predict/{wc.token}/batch_{b}.{wc.vt}.tsv"
                          for b in range(N_BATCHES)],
        manifest=f"{OUT}/predict/_manifests/{{token}}.{{vt}}.gather.txt",
    output:
        table=f"{OUT}/predict/{{token}}.{{vt}}.tsv",
    conda:
        "../envs/traitar3.yaml"
    resources:
        mem_mb=res("gather", "mem_mb"),
        runtime=res("gather", "runtime"),
    shell:
        "python workflow/scripts/gather_predictions.py "
        "--inputs-list {input.manifest} --out {output.table}"


if HAS_SECONDARY:

    rule merge_models:
        """Combine primary (phypat) and secondary (phypat+PGL) into 0/1/2/3 calls."""
        input:
            pmaj=f"{OUT}/predict/{PRIMARY['token']}.majority_vote.tsv",
            psv=f"{OUT}/predict/{PRIMARY['token']}.single_votes.tsv",
            smaj=f"{OUT}/predict/{SECONDARY['token']}.majority_vote.tsv",
            ssv=f"{OUT}/predict/{SECONDARY['token']}.single_votes.tsv",
        output:
            maj=f"{OUT}/predictions_majority-vote_combined.tsv",
            sv=f"{OUT}/predictions_single-votes_combined.tsv",
            fmaj=f"{OUT}/predictions_flat_majority-votes_combined.tsv",
            fsv=f"{OUT}/predictions_flat_single-votes_combined.tsv",
        conda:
            "../envs/traitar3.yaml"
        resources:
            mem_mb=res("merge", "mem_mb"),
            runtime=res("merge", "runtime"),
        params:
            pname=PRIMARY["name"],
            sname=SECONDARY["name"],
            out_dir=lambda wc, output: os.path.dirname(output.maj),
            flat_flag="" if WRITE_FLAT else "--no-flat",
        shell:
            "python workflow/scripts/merge_models.py "
            "--primary-majority {input.pmaj} --primary-single {input.psv} "
            "--secondary-majority {input.smaj} --secondary-single {input.ssv} "
            "--primary-name '{params.pname}' --secondary-name '{params.sname}' "
            "--out-dir {params.out_dir} {params.flat_flag}"

    rule derive_calls:
        """Emit a binary call set (e.g. the phypat+PGL 'reliable' set) from the
        combined 0/1/2/3 matrix. `name` and its codes come from config:derived_call_sets."""
        input:
            combined=f"{OUT}/predictions_majority-vote_combined.tsv",
        output:
            matrix=f"{OUT}/predictions_{{name}}.tsv",
            long=f"{OUT}/predictions_{{name}}_flat.tsv",
        conda:
            "../envs/traitar3.yaml"
        params:
            codes=lambda wc: ",".join(str(c) for c in DERIVED[wc.name]),
            name=lambda wc: wc.name,
        shell:
            "python workflow/scripts/derive_calls.py --combined {input.combined} "
            "--codes {params.codes} --out {output.matrix} --out-long {output.long} "
            "--name {params.name}"
