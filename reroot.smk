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
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium"
    shell:
        """
        python3 tree_time_updated.py -t {input.virus_dir}
        """
