#!/bin/bash

# Set the directory to the fixed path
DIR="../logs/mickg"

# Initialize variables to track the file with the largest number
max_num=-1
max_file=""
max_size=0

# Loop through all files that match the pattern "LNW_cuda-*"
for file in "$DIR"/LNW_cuda-*; do
    # Check if the file matches the pattern "LNW_cuda-*"
    if [[ "$file" =~ LNW_cuda-([0-9]+) ]]; then
        # Extract the number at the end of the filename
        num="${BASH_REMATCH[1]}"

        # Get the size of the file
        size=$(stat --format=%s "$file")

        # Check if this file has the largest number so far
        if (( num > max_num )); then
            max_num=$num
            max_file="$file"
            max_size=$size
        fi
    fi
done

# Output the results
if [[ -n "$max_file" ]]; then
    echo "File with largest number: $max_file"
    echo "Largest number: $max_num"
    echo "File size: $max_size bytes"
else
    echo "No files found matching the pattern 'LNW_cuda-*'."
fi