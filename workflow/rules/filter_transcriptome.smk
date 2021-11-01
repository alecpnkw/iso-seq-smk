rule update_sqlite:
    input:
        "results/sqanti3/{dataset}/{dataset}_classification.txt",
        "results/talon/{dataset}_gene_models.db"
    output:
        touch("results/talon/{dataset}_sqanti3_annotations.done")
    resources:
        mem_mb = 3000
    run:
        import pandas as pd
        import sqlite3

        # connect to SQLite
        con = sqlite3.connect(input[1])
        transcripts = pd.read_sql_query("SELECT transcript_ID, gene_ID from transcripts", con)

        cols = ["isoform", "associated_gene", "associated_transcript", "structural_category", "subcategory", "min_cov"]
        sqanti = pd.read_table(input[0], sep = '\t').loc[cols]
        sqanti_transcripts = pd.merge(transcripts, sqanti, left_on = "transcript_ID", right_on = "isoform", kind = "left")
        sqanti_transcripts.to_sql("sqanti_annotations", con, index = False)

rule filter_transcriptome:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/sqanti3/{dataset}/{dataset}_corrected.gtf",
        "results/talon/{dataset}_sqanti3_annotations.done"
    output:
        "results/filter_transcriptome/{dataset}_whitelist.csv"
    params:
        minFracA = 0.5,
        minCov = 3
    resources:
        mem_mb = 3000
    run:
        import sqlite3
        # using named style
        query = """
        SELECT `gene_ID`, `transcript_ID` FROM transcripts
        WHERE `transcript_ID` NOT IN (
            SELECT DISTINCT `ID`
                FROM `transcript_annotations`
                WHERE (
                    (`attribute` = 'ISM_transcript' AND `value` = 'TRUE') 
                    OR (`attribute` = 'genomic_transcript' AND `value` = 'TRUE')
                )
        ) AND `transcript_ID` IN (
            SELECT DISTINCT `transcript_ID`
                FROM (SELECT `transcript_ID`, COUNT(*) AS `n`
                FROM `observed`
                WHERE (`fraction_As` >= :minFracA)
                GROUP BY `transcript_ID`)
                WHERE (`n` >= :minCov)
        ) AND `transcript_ID` IN (
            SELECT DISTINCT `transcript_ID`
                FROM sqanti_annotations
                WHERE (`min_cov` >= :minCoV)
        )
        """

        # connect to SQLite
        con = sqlite3.connect(input[0])
        cur = con.cursor()

        # execute and write
        with open(output[0], 'w') as o:
            for row in cur.execute(query, {"minFracA": params["minFracA"], "minCov": params["minCov"]}):
                o.write("{0},{1}\n".format(row[0], row[1]))
        
        con.close()

#tabulate
rule filtered_abundance:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/filter_transcriptome/{dataset}_whitelist.csv"
    output: 
        "results/filter_transcriptome/{dataset}_talon_abundance_filtered.tsv" # This name suffix can't be changed!
    conda:
        "../envs/talon.yaml"
    threads:
        4
    resources:
        mem_mb = 5000
    params:
        annot = basename(config["annotation"]),
        output_prefix = lambda wc, output: output[0].split("_abundance_filtered.tsv")[0]
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
rule filtered_GTF:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/filter_transcriptome/{dataset}_whitelist.csv"
    output:
        "results/filter_transcriptome/{dataset}_talon.gtf"
    conda:
        "../envs/talon.yaml"
    params:
        annot = basename(config["annotation"]),
        output_prefix = lambda wc, output: output[0].split("_talon.gtf")[0]
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

rule summarize_transcriptome:
    input:
        "results/filter_transcriptome/{dataset}_talon_abundance_filtered.tsv"
    output:
        "results/filter_transcriptome/{dataset}_summary.nb.html",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/summarize_transcriptome.Rmd"