import re

import sys

# BROUGHT TO US COURTESY OF DEEPSEEK FOR REORDERING NODETERMINISTIC MULTITHREADED OUTPUTS

def reorder_output_file(input_file, output_file):
    # Read all lines from the input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Parse the header, pairs, and footer (elapsed time)
    header = []
    pair_data = []
    footer = []
    current_pair = None
    current_data = []
    
    for line in lines:
        # Check for elapsed time or cleaning up lines
        if line.startswith("Elapsed time") or line.startswith("Cleaning up"):
            if current_pair is not None:
                pair_data.append((current_pair, current_data))
                current_pair = None
                current_data = []
            footer.append(line)
            continue
        
        # Check if this is a pair line (e.g., "3 | 5")
        pair_match = re.match(r'^(\d+) \| (-?\d+)', line)
        if pair_match:
            if current_pair is not None:
                pair_data.append((current_pair, current_data))
            current_pair = int(pair_match.group(1))
            current_data = [line]
        else:
            if line.strip() == '':
                if current_data:  # preserve empty lines within pair data
                    current_data.append(line)
                continue
            if current_pair is None:
                header.append(line)
            else:
                current_data.append(line)
    
    # Add the last pair if exists
    if current_pair is not None:
        pair_data.append((current_pair, current_data))
    
    # Sort the pairs by their number
    pair_data.sort(key=lambda x: x[0])
    
    # Write the output file
    with open(output_file, 'w') as f:
        # Write header
        f.writelines(header)
        
        # Write sorted pairs
        for pair_num, data in pair_data:
            f.writelines(data)
        
        # Write footer (elapsed time and cleaning up)
        f.writelines(footer)
    
    print(f"Output reordered and saved to {output_file}")

# Example usage:
input_file = sys.argv[1]  # Replace with your input file path
output_file = sys.argv[2]   # Replace with desired output file path
reorder_output_file(input_file, output_file)