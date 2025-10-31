#!/bin/bash

for dir in test_CA*/; do
    # Check if it's actually a directory
    if [ -d "$dir" ]; then
        # Check if any .mp4 files exist inside the directory
        if ! compgen -G "$dir*.mp4" > /dev/null; then
            echo "No .mp4 found in $dir, running Analysis.py"
            # Run Analysis.py in that directory
            (cd "$dir" && python3 Analysis.py)
        else
            echo ".mp4 file exists in $dir, skipping"
        fi
    fi
done
