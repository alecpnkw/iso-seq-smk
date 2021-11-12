# conditional rule definition based on config...
if config.get("short_read_splice_junctions", None) is not None:
    rule merge_sjs:
        input:  
            config["short_read_splice_junctions"]
        output:
            "results/merge_sjs/{dataset}/SJ.merged.tab"
        conda: 
            "../envs/R.yaml"
        script:
            "../scripts/merge_sjs.R"


# this looks like it pulls the first field of the chr names...
rule rename_fasta_chromosomes:
    input: config["reference"]
    output: "results/renamed_ref/{reference}.fa"
    shell:
        """
        cut -d ' ' -f1 {input} > {output}
        """

# start of TranscriptClean
# one .sam at a time
def transcriptclean_cmd(config):
    cmd_str = """
        python3  {params.tc_exec} \
        --sam {input[0]} \
        --genome {input[1]} \
        --tmpDir {params.tempdir} \
        --threads {threads} \
        --deleteTmp \
        --canonOnly \
        --primaryOnly \
        --outprefix {params.output_prefix}"""
    if config.get("short_read_splice_junctions", None) is not None:
        cmd_str += " --spliceJns {input[2]}"
    return cmd_str

def TC_INPUT(wc):
    inputs = [
        "results/minimap2_filt/{sample}.sam",
        "results/renamed_ref/{0}".format(basename(config["reference"]))
        ]
    if config.get("short_read_splice_junctions", None) is not None:
        inputs.append("results/merge_sjs/{dataset}/SJ.merged.tab")
    return inputs

TC_SHELL = transcriptclean_cmd(config)
rule transcript_clean:
    input:
        TC_INPUT
    output:
        multiext("results/transcript_clean/{dataset}/{sample}/TC_clean", ".sam", ".log", ".fa", ".TE.log")
    params:
        output_prefix = lambda wc, output: dirname(output[0]),
        tc_exec = config["transcriptclean_path"],
        tempdir = "tmp/{dataset}_{sample}/"
    threads:
        10
    resources:
        mem_mb = 8000
    conda: 
        "../envs/transcript_clean.yaml"
    shell:
        TC_SHELL