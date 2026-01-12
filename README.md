This pipeline is intended to run treetime on viral usher trees. although tree time results are mixe

# Set up 
treetime is not run directly, instead it is run through a script called run_treetime.py developed by Angie Hinrichs and modified slightly by Lily Karim for automation 

Before running run_treetime.py, viral_usher_trees (put link here) should be cloned into the main dir called reroot_pipeline. this clone can and should be maintained each month through a git pull which will retrieve the most up to date trees. this repo should not be added to any git pushes. after viral_usher_trees is cloned into reroot_pipeline, the script get_fastas.sh should be moved into viral_usher_trees and run to retrieve fasta info for run_treetime. 

# To get tree time to run please create a conda environment from the yml file in envs.  


# For a single tree, tree time can be run with the following command 
`python3 run_treetime.py -t {path to virus directory}` note that the script does not assume your directory but it will assume that all output is written into the viral directory 