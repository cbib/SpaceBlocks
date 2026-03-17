rule spaceranger_count:
    """
    Run Space Ranger ``count`` on a single Visium HD sample.

    Space Ranger must be installed locally — its licence does not permit
    redistribution in container images.  The executable path is set in
    ``config.yaml`` under the ``spaceranger`` key.

    Real output files (web_summary.html and the raw bin matrix) are NOT tracked
    as the output of this rule due to Martian restrictions. These files are
    instead tracked in downstream rules, giving Snakemake proper provenance.
    """
    input:
        fastqs=lambda wc: f"{SAMPLES[wc.sample]['fastq_dir']}/fastqs/",
        cytaimage=lambda wc: (
            f"{SAMPLES[wc.sample]['fastq_dir']}/images/{wc.sample}_cyta.tiff"
        ),
        image=lambda wc: (
            f"{SAMPLES[wc.sample]['fastq_dir']}/images/{wc.sample}.tiff"
        ),
        loupe_alignment=lambda wc: (
            f"{SAMPLES[wc.sample]['fastq_dir']}/images/{wc.sample}.json"
        ),
        transcriptome=config["transcriptome"],
        probe_set=config["probe_set"],
    output:
        done_flag=f"{OUTDIR_SR}/{{sample}}/.done",
    params:
        outdir=f"{OUTDIR_SR}/{{sample}}",
        slide=lambda wc: SAMPLES[wc.sample]["slide"],
        area=lambda wc: SAMPLES[wc.sample]["area"],
        spaceranger=config["spaceranger"],
        fastqs_formatted=fastq_dirs_comma_separated,
    log:
        out=f"{LOGDIR}/spaceranger_count/{{sample}}.out",
        err=f"{LOGDIR}/spaceranger_count/{{sample}}.err",
    benchmark:
        f"{LOGDIR}/benchmarks/spaceranger_count/{{sample}}.tsv"
    threads:
        get_resource("spaceranger_count", "threads")
    resources:
        mem_mb=get_resource("spaceranger_count", "mem_mb"),
        runtime=get_resource("spaceranger_count", "runtime"),
    shell:
        """
        (
        echo "=== Space Ranger version ==="
        {params.spaceranger} --version

        {params.spaceranger} count \
            --id={wildcards.sample} \
            --transcriptome={input.transcriptome} \
            --fastqs={params.fastqs_formatted} \
            --sample={wildcards.sample} \
            --cytaimage={input.cytaimage} \
            --image={input.image} \
            --slide={params.slide} \
            --area={params.area} \
            --loupe-alignment={input.loupe_alignment} \
            --probe-set={input.probe_set} \
            --output-dir={params.outdir} \
            --create-bam=false
        ) > {log.out} 2> {log.err}

        touch {output.done_flag}
        """
