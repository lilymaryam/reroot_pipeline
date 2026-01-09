This pipeline is intended to run treetime on viral usher trees. although tree time results are mixe

# Set up 
treetime is not run directly, instead it is run through a script called tree_time_updated.py developed by Angie Hinrichs and modified slightly by lily karim for automation 

Before running tree_time_updated, viral_usher_trees (put link here) should be cloned into the main dir called reroot_pipeline. this clone can and should be maintained each month through a git pull whcich will retrieve the most up to date trees. this repo should not be added to any git pushes. after viral_usher_trees is cloned into reroot_pipeline, the script get_fastas.sh should be moved into viral_usher_trees and run to retrieve fasta info for tree_time_updated. 

# To get tree time to run please create a conda environment from the yml file in envs.  


# For a single tree, tree time can be run with the following command 
`python3 tree_time_updated.py -t {path to virus directory}` note that the script does not assume your directory but it will assume that all output is written into the viral directory 