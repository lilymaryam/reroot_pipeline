import os

#make this a variable?
tree_dir = "../viral_usher_trees/trees/"

# Only keep directories that contain a .pb.gz file with the same name

trees = []
for d in os.listdir(tree_dir):
    trees.append(d)

print(trees)
 

rule all:
    input:
        expand("outputs/{tree}.facts", tree=trees)

rule get_trees:
    input:
        f"{tree_dir}" + "{tree}/optimized.pb.gz"
    output:
        "outputs/{tree}.facts"
    shell:
        "matUtils summary -i {input} > {output}"

rule make_outputs_dir:
    output:
        directory("outputs")
    shell:
        "mkdir -p outputs"

rule reroot:
    input:
        tree=f"{tree_dir}" + "{tree}/optimized.pb.gz",
        metadata=f"{tree_dir}" + "{tree}/metadata.tsv.gz"
    output:
