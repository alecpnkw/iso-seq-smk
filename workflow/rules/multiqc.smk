rule multiqc:
    input:
        "results/pbmm2/mapped.bam"
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