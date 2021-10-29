rule filter_transcriptome:
    input:
        "results/talon/{dataset}_gene_models.db",
        "results/sqanti3/{dataset}/{dataset}_corrected.gtf"
    output:
        "results/filter_transcriptome/{dataset}_whitelist.csv"
    params:
        minFracA = 0.5,
        minCov = 3
    run:
        import sqlite3
        # using qmark (?) style
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
                WHERE (`fraction_As` >= ?)
                GROUP BY `transcript_ID`)
                WHERE (`n` >= ?)
        )"""

        # connect to SQLite
        con = sqlite3.connect(input["db"])
        cur = con.cursor()

        # execute and write
        with open(snakemake.output[0], 'w') as o:
            for row in cur.execute(query, (params["minFracA"], params["min_cov"])):
                o.write(",".join(row))