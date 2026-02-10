#!/bin/bash

# Define the base source and destination paths
SRC_DIR="viral_usher_trees/trees"
DEST_BASE="treetime_feb"

# Create the main destination directory if it doesn't exist
mkdir -p "$DEST_BASE"

# Iterate through every directory in the source path
for dir_path in "$SRC_DIR"/*; do
    # Ensure we are working with a directory
    if [ -d "$dir_path" ]; then
        # Get the directory name without the full path
        dir_name=$(basename "$dir_path")
        
        # Create the corresponding directory in treetime_feb
        new_dest="$DEST_BASE/$dir_name"
        mkdir -p "$new_dest"

        # 1. Copy treetime.log if it exists
        if [ -f "$dir_path/treetime.log" ]; then
            cp "$dir_path/treetime.log" "$new_dest/"
        fi

        # 2. If treetime_out exists, copy the whole directory
        if [ -d "$dir_path/treetime_out" ]; then
            cp -r "$dir_path/treetime_out" "$new_dest/"
            
            # 3. If treetime_out exists, also copy specific files to the {dir_name} root
            # Using find or glob to match the files
            cp "$dir_path"/timetree_rerooted.*gz "$new_dest/" 2>/dev/null
            cp "$dir_path"/treetime_rerooted_*.fa "$new_dest/" 2>/dev/null
        fi
    fi
done

echo "Transfer complete."
