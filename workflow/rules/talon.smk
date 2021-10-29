# initialize db (uses SQLite)...
# follow up with Denis about additional options...
rule talon_initialize:
    input:
        config["annotation"] #pass in via config
    output:
        "results/talon/{dataset}_gene_models.db", #only one genemodels.db per run...
    conda:
        "../envs/talon.yaml"
    params:
        annot = basename(config["annotation"]),
        output_prefix = lambda wildcards, output: output[0].split('.db')[0]
    threads:
        4
    resources:
        mem_mb = 6000 # change to mem_mb
    shell: 
        """
        talon_initialize_database \
        --f {input} \
        --a {params.annot} \
        --g {wildcards.dataset} \
        --o {params.output_prefix}
        """

#label/flag reads...
rule talon_label:
    input:
        "results/transcript_clean/{dataset}/{sample}/TC_clean.sam" #aligned, filtered reads
    output:
        "results/talon_label/{dataset}/{sample}_labeled.sam",
        "results/talon_label/{dataset}/{sample}_read_labels.tsv"
    conda:
        "../envs/talon.yaml"
    threads: 2
    params:
        output_prefix = lambda wc, output: output[0].split('_labeled.sam')[0],
        ref = config["reference"],
        ar = 20,
        tempdir = "tmp/talon/{sample}/"
    shell:
        """
        talon_label_reads \
        --f {input} \
        --g {params.ref}  \
        --t {threads} \
        --ar {params.ar} \
        --tmpDir {params.tempdir} \
        --deleteTmp \
        --o {params.output_prefix}
        """

#writes a new config.csv... 
#tweak this... 
rule create_talon_config:
    input:
        # all labeled files as input... 
        expand("results/talon_label/{dataset}/{sample}_labeled.sam", zip, sample = samples.index, dataset = samples.dataset)
    params:
        config = config["samples"]
    output:
        "results/talon/{dataset}_config.csv"
    run:
        import pandas as pd
        tbl = pd.read_csv(params.config, sep = '\t')
        tbl["file"] = input
        keeps = tbl.dataset == wildcards.dataset
        cols = ["sample", "dataset", "source", "file"]
        tbl.loc[keeps, cols].to_csv(output[0], header = False, index = False)

#TALON... 
rule talon:
    input:
        config = "results/talon/{dataset}_config.csv",
        db = "results/talon/{dataset}_gene_models.db"
    params: 
        build = basename(config["reference"])
    threads:
        10
    resources:
        mem_mb = 10000
    output:
        "results/talon/{dataset}_QC.log"
    conda:
        "../envs/talon.yaml"
    shell: 
        """
        talon \
        --f {input.config} \
        --db {input.db} \
        --build {wildcards.dataset} \
        --threads {threads} \
        --cov 0.99 \
        --identity 0.95 \
        --o results/talon/{wildcards.dataset}
        """

# transcript model filtering... 
# dataset names need to match talon config...
# this differs from Denis's run... which supplies an SQL query
rule talon_filter:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/talon/{dataset}_QC.log"
    output:
        "results/talon/{dataset}_filtered_transcripts.csv"
    params:
        samples = lambda wc: ",".join(samples.index[samples.dataset == wc.dataset]),
        maxFracA = 0.5,
        minCount = 0,
        minDatasets = 0,
        annot = basename(config["annotation"])
    threads:
        4
    resources:
        mem_mb = 5000 
    conda:
        "../envs/talon.yaml"
    shell:
        """
        talon_filter_transcripts \
        --db {input[0]} \
        --datasets {params.samples} \
        -a {params.annot} \
        --maxFracA {params.maxFracA} \
        --minCount {params.minCount} \
        --minDatasets {params.minDatasets} \
        --o {output}
        """

#tabulate
rule talon_abundance:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/talon/{dataset}_filtered_transcripts.csv"
    output: 
        "results/talon/{dataset}_talon_abundance_filtered.tsv"
    conda:
        "../envs/talon.yaml"
    threads:
        4
    resources:
        mem_mb = 5000
    params:
        annot = basename(config["annotation"]),
        output_prefix = lambda wc, output: output[0].split("_talon_abundance_filtered.tsv")[0]
    shell:
        """
        talon_abundance \
        --db {input[0]} \
        --whitelist {input[1]} \
        --annot {params.annot} \
        --build {wildcards.dataset} \
        --o {params.output_prefix}
        """

#GTF
rule talon_GTF:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/talon/{dataset}_filtered_transcripts.csv"
    output:
        "results/talon/{dataset}_talon.gtf"
    conda:
        "../envs/talon.yaml"
    params:
        annot = basename(config["annotation"]),
        output_prefix = lambda wc, output: output[0][:-10]
    threads:
        4
    resources:
        mem_mb = 5000
    shell: 
        """
        talon_create_GTF \
        --db {input[0]} \
        --whitelist {input[1]} \
        -a {params.annot} \
        --build {wildcards.dataset} \
        --o {params.output_prefix}
        """