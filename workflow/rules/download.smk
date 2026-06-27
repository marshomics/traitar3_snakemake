# Pfam-A database acquisition and (optional) subsetting to the model feature set.

rule download_pfam:
    """Download + gunzip Pfam-A 27.0. `protected` guards the ~1.2 GB result."""
    output:
        hmm=protected(config["pfam"]["hmm"]),
    params:
        url=config["pfam"]["url"],
    conda:
        "../envs/traitar3.yaml"
    resources:
        mem_mb=res("download", "mem_mb"),
        runtime=res("download", "runtime"),
    shell:
        r"""
        mkdir -p "$(dirname {output.hmm})"
        curl -fL --retry 8 --retry-delay 10 -o {output.hmm}.gz "{params.url}"
        gunzip -c {output.hmm}.gz > {output.hmm}
        rm -f {output.hmm}.gz
        # sanity check: must contain HMMER profiles
        grep -qm1 '^HMMER3' {output.hmm} || (echo "ERROR: downloaded file is not an HMMER db" >&2; exit 1)
        """


rule model_pfam_accessions:
    """Union of Pfam accessions referenced by the configured phenotype models."""
    input:
        archives=[m["archive"] for m in MODELS],
    output:
        accs=f"{OUT}/db/model_pfam_accessions.txt",
    conda:
        "../envs/traitar3.yaml"
    shell:
        "python workflow/scripts/model_pfam_accessions.py {input.archives} -o {output.accs}"


rule subset_pfam:
    """Reduce Pfam-A to just the families the models use (identical hits, faster)."""
    input:
        accs=f"{OUT}/db/model_pfam_accessions.txt",
        hmm=config["pfam"]["hmm"],
    output:
        hmm=config["pfam"]["subset_hmm"],
    conda:
        "../envs/traitar3.yaml"
    resources:
        mem_mb=res("subset_db", "mem_mb"),
        runtime=res("subset_db", "runtime"),
    shell:
        "python workflow/scripts/subset_hmm.py {input.accs} {input.hmm} {output.hmm}"
