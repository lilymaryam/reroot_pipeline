configfile: "config.yaml"

trees = config["viruses"]
#note this must be the path all the way to viral_usher_trees/trees
path_to_viruses = config["path_to_viruses"]

rule all:
    input:
        expand(f"{path_to_viruses}/{{tree}}/treetime.log", tree=trees)

rule reroot:
    input:
        virus_dir=f"{path_to_viruses}/{{tree}}/"
    output:
        treetime_log=f"{path_to_viruses}/{{tree}}/treetime.log"
    conda: "reroot_env"
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium",
        slurm_extra="--export=ALL --mail-type=FAIL --mail-user=lkarim@ucsc.edu",
    shell:
        """
        python3 run_treetime.py -t {input.virus_dir}
        """
