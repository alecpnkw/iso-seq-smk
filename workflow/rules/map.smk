# FLNC .bam to .fastq
rule flnc_to_fastq:
    input: INPUT_BAM
    output: "results/flnc_fastq/{sample}.fastq.gz"
    #conda:
    #    "../envs/samtools.yaml"
    envmodules:
        "samtools/1.11"
    shell:
        """
        samtools view {input} | \
        awk '{{printf("@%s\\n%s\\n+\\n%s\\n", $1, $10, $11)}}' | \
        gzip > {output}
        """
        
# Todo: add minimap2 conda env
rule minimap2_sam:
    input:
        target=config["reference"], # can be either genome index or genome fasta
        query="results/flnc_fastq/{sample}.fastq.gz"
    output:
        "results/minimap2/{sample}.sam" # can add additional output files
    params:
        output_prefix = lambda wc, output: output[0].split(".sam")[0],
        threads = 6
    threads: 6
    resources:
        mem_mb = 6000 #change to mem_mb
    envmodules:
        "samtools/1.9",
        "minimap2/2.17"
    shell:
        """
        minimap2 -ax splice -uf --secondary=no -C5 -t {params.threads} --MD {input.target} {input.query} | \
        samtools sort -O sam -o {output} && \
        samtools flagstat {output} > {params.output_prefix}.flagstat && \
        samtools idxstats {output} > {params.output_prefix}.idxstats && \
        samtools stats {output} > {params.output_prefix}.stats && \
        samtools view {output} | cut -f1,5 > {params.output_prefix}.mapq
        """

rule filter_sam:
    input:
        "results/minimap2/{sample}.sam"
    output:
        "results/minimap2_filt/{sample}.sam"
    threads:
        5
    params:
        output_prefix = lambda wc, output: output[0].split('.sam')[0],
        threads = 30
    #conda:
    #    "../envs/sambamba.yaml"
    envmodules:
        "samtools/1.9",
        "sambamba/0.5.6"
    shell:
        """
        sambamba view \
        --with-header \
        --nthreads {params.threads} \
        --sam-input \
        --format sam \
        --filter "not unmapped and mapping_quality >= 50" \
        {input} | \
        samtools sort -O sam -o {output} && \
        samtools flagstat {output} > {params.output_prefix}.flagstat && \
        samtools idxstats {output} > {params.output_prefix}.idxstats && \
        samtools stats {output} > {params.output_prefix}.stats && \
        samtools view {output} | cut -f1,5 > {params.output_prefix}.mapq
        """