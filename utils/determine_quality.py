"""
Title: determine_quality.py
Function: This script is used to determine the quality of models based on sampling results.
Date: 11-11-2024
Author: Nils Smit
"""

from sys import argv

def parse_sampling_results(sampling_file):
    """Parse the sampling results file and store the results in a dictionary
    
    Args:
        sampling_file (str): Path to the sampling results file
    
    Returns:
        dict: Dictionary containing the results for each model
    """
    with open(sampling_file) as f:
        lines = f.readlines()

    # Initialize the result dictionary
    res_dict = {}

    # Process each line in the file
    for line in lines:
        # Split line by tabs
        splits = line.strip().split("\t")
        id = splits[0]
        
        # Convert each tuple string to a tuple and reformat it
        tuples = [eval(t) for t in splits[1:]]
        results = []
        for tup in tuples:
            (lrmsd, irmsd, fnat) = tup
            results.append((fnat, irmsd, lrmsd))  # Store as (fnat, irmsd, lrmsd)
        
        # Add to dictionary
        res_dict[id] = results

    return res_dict

# Quality check functions
def check_fnat(fnat):
    """Check the FNAT value and return a score based on the following criteria:

    Args:
        fnat (float): FNAT value to check

    Returns:
        int: Score based on FNAT value
    """
    if fnat >= 0.5:
        return 3
    elif fnat >= 0.3:
        return 2
    elif fnat >= 0.1:
        return 1
    else:
        return 0

def check_irmsd(irmsd):
    """Check the IRMSD value and return a score based on the following criteria:
    
    Args:
        irmsd (float): IRMSD value to check
    
    Returns:
        int: Score based on IRMSD value
    """
    if irmsd <= 1.0:
        return 3
    elif irmsd <= 2.0:
        return 2
    elif irmsd <= 4.0:
        return 1
    else:
        return 0

def check_lrmsd(lrmsd):
    """Check the LRMSD value and return a score based on the following criteria:
    
    Args:
        lrmsd (float): LRMSD value to check
        
    Returns:
        int: Score based on LRMSD value
    """
    if lrmsd <= 1.0:
        return 3
    elif lrmsd <= 5.0:
        return 2
    elif lrmsd <= 10.0:
        return 1
    else:
        return 0

def determine_model_quality(value):
    """Determine the overall quality of a model based on FNAT, IRMSD, and LRMSD scores

    Args:
        value (tuple): Tuple containing FNAT, IRMSD, and LRMSD scores
    
    Returns:
        str: Quality level of the model (high, medium, acceptable, incorrect)
    """
    (fnat, irmsd, lrmsd) = value
    fnat_score = check_fnat(fnat)
    irmsd_score = check_irmsd(irmsd)
    lrmsd_score = check_lrmsd(lrmsd)
    
    if fnat_score == 3 and (irmsd_score == 3 or lrmsd_score == 3):
        return 'high'
    elif fnat_score >= 2 and (irmsd_score >= 2 or lrmsd_score >= 2):
        return 'medium'
    elif fnat_score >= 1 and (irmsd_score >= 1 or lrmsd_score >= 1):
        return 'acceptable'
    else:
        return 'incorrect'

def write_results_to_tsv(res_dict, output_file):
    """Write model qualities in a tabular format with ranks as rows and models as columns

    Args:
        res_dict (dict): Dictionary with model names as keys and quality data as values
        output_file (str): Path to output TSV file
    """
    # Find the maximum number of ranks for formatting
    max_ranks = max(len(qualities) for qualities in res_dict.values())

    # Open file to write results
    with open(output_file, 'w') as f:
        # Write the header with each model key followed by "_quality"
        header = "\t".join([model_key for model_key in res_dict.keys()])
        f.write(f"\t{header}\n")  # Initial tab for aligning rank

        # Write each rank row
        for rank in range(max_ranks):
            # Start each line with "rank_X" (e.g., rank_1, rank_2, ...)
            row = [f"rank_{rank + 1}"]

            # Append quality for each model at the current rank, or empty if it doesn't exist
            for model_key in res_dict.keys():
                if rank < len(res_dict[model_key]):
                    quality = determine_model_quality(res_dict[model_key][rank])
                else:
                    quality = ""  # No value for this rank if model has fewer entries
                row.append(quality)

            # Write the row to the file
            f.write("\t".join(row) + "\n")
    
if __name__ == "__main__":
    input_file = argv[1]
    output_file = argv[2]  # Define output TSV file path
    
    res_dict = parse_sampling_results(input_file)
    
    # Write to TSV
    write_results_to_tsv(res_dict, output_file)

    print(f"Results written to {output_file}")