rule update_sqlite:
    input:
        "results/sqanti3/{dataset}/{dataset}_classification.txt",
        "results/talon/{dataset}_gene_models.db"
    output:
        touch("results/talon/{dataset}_sqanti3_annotations.done")
    resources:
        mem_mb = 5000
    run:
        import pandas as pd
        import sqlite3

        # connect to SQLite
        con = sqlite3.connect(input[1])
        transcripts = pd.read_sql_query("SELECT `ID`, `value` AS `transcript_ID` FROM `transcript_annotations` WHERE (`attribute` = 'transcript_id')", con)

        cols = ["isoform", "associated_gene", "associated_transcript", "structural_category", "subcategory", "min_cov"]
        sqanti = pd.read_table(input[0], sep = '\t').loc[:,cols]
        sqanti_transcripts = pd.merge(transcripts, sqanti, left_on = "transcript_ID", right_on = "isoform", how = "left")
        sqanti_transcripts.to_sql("sqanti_annotations", con, index = False)
        con.close()

rule filter_transcriptome:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/sqanti3/{dataset}/{dataset}_corrected.gtf",
        "results/talon/{dataset}_sqanti3_annotations.done"
    output:
        "results/filter_transcriptome/{dataset}_whitelist.csv"
    resources:
        mem_mb = 3000
    run:
        import sqlite3
        # using named style
        query = """
        SELECT gene_ID, t.transcript_ID, transcript_novelty, n_exons, min_cov, n, avg_fraction_As
        FROM transcripts t
        LEFT JOIN (SELECT ID, CASE 
        WHEN sum(`attribute` = 'transcript_status' AND `value` = 'KNOWN') > 0 THEN 'KNOWN'
        WHEN sum(`attribute` = 'ISM_transcript' AND `value` = 'TRUE') > 0 THEN 'ISM'
        WHEN sum(`attribute` = 'NIC_transcript' AND `value` = 'TRUE') > 0 THEN 'NIC'
        WHEN sum(`attribute` = 'NNC_transcript' AND `value` = 'TRUE') > 0 THEN 'NNC'
        WHEN sum(`attribute` = 'intergenic_transcript' AND `value` = 'TRUE') > 0 THEN 'intergenic'
        WHEN sum(`attribute` = 'antisense_transcript' AND `value` = 'TRUE') > 0 THEN 'antisense'
        WHEN sum(`attribute` = 'genomic_transcript' AND `value` = 'TRUE') > 0 THEN 'genomic'
        ELSE 'OTHER_NOVEL' END 
        AS `transcript_novelty`
        FROM transcript_annotations
        GROUP BY `ID`) a ON t.transcript_ID=a.ID
        LEFT JOIN sqanti_annotations s ON t.transcript_ID=s.ID
        LEFT JOIN (SELECT t.transcript_ID, COUNT(*) as `n`, AVG(o.fraction_As) as avg_fraction_As
            FROM transcripts t
            LEFT JOIN observed o ON o.transcript_ID=t.transcript_ID
            WHERE o.fraction_As <= 0.5
            GROUP BY t.transcript_ID, `n_exons`) o ON t.transcript_ID=o.transcript_ID
        WHERE `transcript_novelty` NOT IN ('genomic', 'ISM') 
        AND (`min_cov` >= 5 OR `n_exons` = 1)
        AND (
            (`n` >= 1 AND `transcript_novelty` == 'KNOWN') 
            OR (`n` >= 2 AND `transcript_novelty` IN ('NIC', 'NNC'))
            OR (`n` >= 3 AND `transcript_novelty` IN ('intergenic', 'antisense'))
            )
        """

        # connect to SQLite
        con = sqlite3.connect(input[0])
        cur = con.cursor()

        # execute and write
        with open(output[0], 'w') as o:
            for row in cur.execute(query):
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