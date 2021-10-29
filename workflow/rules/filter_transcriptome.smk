rule filter_transcriptome:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/sqanti3/{dataset}/{dataset}_corrected.gtf"
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
        )"""

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
        "results/filter_transcriptome/{dataset}_abundance.tsv"
    conda:
        "../envs/talon.yaml"
    threads:
        4
    resources:
        mem_mb = 5000
    params:
        annot = basename(config["annotation"]),
        output_prefix = lambda wc, output: output[0].split("_abundance.tsv")[0]
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