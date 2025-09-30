import os

#make this a variable?
tree_dir = "viral_usher_trees/trees/"

# Only keep directories that contain a .pb.gz file with the same name

trees = []
for d in os.listdir(tree_dir):
    trees.append(d)

rule all:
    input:
        expand("viral_usher_trees/trees/{tree}/treetime_out/rerooted.newick", tree=trees)

rule reroot:
    output:
        "viral_usher_trees/trees/{tree}/treetime_out/rerooted.newick"
    params:
        tree_dir = tree_dir
    resources:
        mem_mb=4000,
        runtime=720,
        slurm_partition="medium"
    shell:
        """
        python3 tree_time.py -v {wildcards.tree} -d {params.tree_dir}/
        """