rule multiqc:
    input:
        expand("results/talon/{dataset}_talon.gtf", dataset = samples["dataset"]),
        expand("results/talon/{dataset}_talon_abundance_filtered.tsv", dataset = samples["dataset"]),
    output:
        "results/multiqc/multiqc-report.html"
    conda: 
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc \
        --no-data-dir \
        --filename {output} \
        results/
        """