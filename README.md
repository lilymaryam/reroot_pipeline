# Set up 
treetime is not run directly, instead it is run through a script called run_treetime.py developed by Angie Hinrichs and modified slightly by Lily Karim for automation 

Before running run_treetime.py, viral_usher_trees (put link here) should be cloned into the main dir called reroot_pipeline. this clone can and should be maintained each month through a git pull which will retrieve the most up to date trees. this repo should not be added to any git pushes. after viral_usher_trees is cloned into reroot_pipeline, the script get_fastas.sh should be moved into viral_usher_trees and run to retrieve fasta info for run_treetime. 

run_treetime.py relies on a script called alter_gbff.py which is also developed by angie hinrichs. I have a version of the script available in reroot_pipeline but this can also be downloaded from viral_usher_trees/scripts

**To get treetime to run please create a conda environment from the yml file in envs.** 

In addition to all dependencies in `envs/environment.yaml` this process also relies on usher/matUtils being available as a command. I do not use the conda version of usher, I have a conda local build of usher that I added to my $PATH. Before using this pipeline make sure you have a version of usher that can be called with `usher` or `matUtils`. This program relies on commit a87e940f of usher. I used the local build (note local build instruction are insufficient, i need to do something about that).  



## For a single tree, tree time can be run with the following command 
`python3 run_treetime.py -t {path to virus directory}` note that the script does not assume your directory but it will assume that all output is written into the viral directory 

## To run treetime on all viral samples use reroot.smk to "parallelize"
To run with slurm on a cluster run the command `snakemake -s reroot.smk --slurm --jobs 20 --use-conda --rerun-incomplete -k`. `--use-conda` may not be necessary as the slurm command will `--export=ALL`. reroot.smk has resources delegated in the rules. these work for my cluster defaults but may need to be modified for use elsewhere. If matUtils is in HOME path export all should allow it to run in cluster. 

To run on a server use the command `snakemake -s reroot.smk --cores 5 --rerun-incomplete -k`. this will run rules based on the number of cores alloted which will also improve the efficiency of processing all of the data. 

Current version of pipeline will put all timetree reroots directly into viral usher trees. this is fine but may be less helpful if running over and over. To move files out of viral usher trees copyteetime.sh can be used to transfer files elsewhere. 