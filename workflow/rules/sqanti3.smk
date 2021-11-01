# Config fields used: 
#   - cDNA_cupcake_path
#   - sqanti_path
#   - annotation
#   - CAGE_peaks
#   - polyA_motifs
#   - short_read_splice_junctions

rule cDNA_cupcake_init:
    output: touch("results/cDNA_cupcake_build.done")
    params: 
        path = config["cDNA_cupcake_path"]
    conda: "../envs/sqanti3.yaml"
    shell: 
        """
        cd {params.path}
        python {params.path}/setup.py build
        python {params.path}/setup.py install
        """

def PARSE_SJ(wc): 
    if config.get("short_read_splice_junctions", None) is not None:
        return ",".join(config["short_read_splice_junctions"])
    else:
        return None

def sqanti_cmd(config):
    cmd_str = """
    export PYTHONPATH={params.cupcake_path}/sequence/
    export PYTHONPATH=$PYTHONPATH:{params.cupcake_path}/
    python {params.sqanti}/sqanti3_qc.py \
    {input.isoforms} {input.annot} {input.genome} \
    --cage_peak {params.cage_path} \
    --polyA_motif_list {params.polya_motifs} \
    -d {params.out_dir} \
    -o {wildcards.dataset} \
    {params.addtl_opt}\
    """
    if config.get("short_read_splice_junctions", None) is not None:
        cmd_str += "-c {params.sj_files} "
    return cmd_str

SQANTI_SHELL = sqanti_cmd(config)

rule sqanti_qc:
    input:
        "results/cDNA_cupcake_build.done",
        isoforms = "results/talon/{dataset}_talon.gtf",
        annot = config["annotation"],
        genome = lambda wc: "results/renamed_ref/{0}".format(basename(config["reference"])) # pass in renamed chr for consistency
    output:
        "results/sqanti3/{dataset}/{dataset}_corrected.gtf",
        "results/sqanti3/{dataset}/{dataset}_classification.txt"
    conda:
        "../envs/sqanti3.yaml"
    params:
        sqanti = config["sqanti_path"],
        addtl_opt = "--force_id_ignore --skipORF --report pdf ",
        sj_files = PARSE_SJ,
        out_dir = lambda wc, output: dirname(output[0]),
        cage_path = config["CAGE_peaks"],
        polya_motifs = config["polyA_motifs"],
        cupcake_path = config["cDNA_cupcake_path"]
    threads: 8
    resources:
        mem_mb = 8000
    shell:
        SQANTI_SHELL