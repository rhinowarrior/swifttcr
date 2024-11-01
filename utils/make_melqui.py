# -*- coding: utf-8 -*-
"""
Name: make_melqui.py
Function: Create a melqui plot from the rmsd values of the models
Date: 04-06-2024 11:47
Author: Yanic
"""

"""
Input:
    - A file containing rmsd values for each model
    - A file containing the output pdf file

Example usage:
python3 utils/make_melqui.py input_rmsd_values.txt output.pdf
"""

#imports
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
import numpy as np
from sys import argv

def check_fnat(fnat):
    """Function to check if fnat is high, medium, acceptable, incorrect coded as 3, 2, 1, 0
    
    Args:
        fnat (float): Fraction of native contacts
    return:
        int: 3 if high, 2 if medium, 1 if acceptable, 0 if incorrect
    """
    if fnat >= 0.5:
        return 3
    elif fnat >= 0.3:
        return 2
    elif fnat >= 0.1:
        return 1
    elif fnat < 0.1:
        return 0
    else:
        raise Exception("Invalid fnat value.")
    
def check_irmsd(irmsd):
    """Function to check if irmsd is high, medium, acceptable, incorrect coded as 3, 2, 1, 0
    
    Args:
        irmsd (float): Interface RMSD
        
    return:
        int: 3 if high, 2 if medium, 1 if acceptable, 0 if incorrect
    """
    if irmsd <= 1.0:
        return 3
    elif irmsd <= 2.0:
        return 2
    elif irmsd <= 4.0:
        return 1
    elif irmsd > 4.0:
        return 0
    else:
        raise Exception("Invalid irmsd value.")

def check_lrmsd(lrmsd):
    """Function to check if lrmsd is high, medium, acceptable, incorrect coded as 3, 2, 1, 0

    Args:
        lrmsd (float): Ligand RMSD
    
    return:
        int: 3 if high, 2 if medium, 1 if acceptable, 0 if incorrect
    """
    if lrmsd <= 1.0:
        return 3
    elif lrmsd <= 5.0:
        return 2
    elif lrmsd <= 10.0:
        return 1
    elif lrmsd > 10.0:
        return 0
    else:
        raise Exception("Invalid lrmsd value.")


def _color_code(value, rigid):
    """Color code the models based on the rmsd values
    
    Args:
        value (tuple): Tuple containing the fnat, irmsd, and lrmsd values
        rigid (bool): True if rigid, False otherwise
        
    return:
        str: color code for the model
    """
    # Unpack the tuple into fnat, irmsd, and lrmsd
    (fnat, irmsd, lrmsd) = value
    fnat_score = check_fnat(fnat)
    irmsd_score = check_irmsd(irmsd)
    lrmsd_score = check_lrmsd(lrmsd)
    if rigid:
        if fnat_score == 3 and (irmsd_score == 3 or lrmsd_score == 3):
            return '#003600'  
        elif fnat_score >= 2 and (irmsd_score >= 2 or lrmsd_score >= 2):
            return 'green'
        elif fnat_score >= 1 and (irmsd_score >= 1 or lrmsd_score >= 1):
            return 'lightgreen'
        else:
            return 'white'
    else:
        if fnat_score == 3 and (irmsd_score == 3 or lrmsd_score == 3):
            return 'midnightblue'
        elif fnat_score >= 2 and (irmsd_score >= 2 or lrmsd_score >= 2):
            return 'blue'
        elif fnat_score >= 1 and (irmsd_score >= 1 or lrmsd_score >= 1):
            return 'lightblue'
        else:
            return 'white'


def parse_sampling_results(sampling_file):
    """Parse the sampling results file and store the results in a dictionary
    
    Args:
        sampling_file (str): Path to the sampling results file
    
    return:
        dict: Dictionary containing the results for each model
    """
    f = open(sampling_file)
    lines = f.readlines()
    f.close()
    res_dict = {}
    for line in lines:
        splits = line.split("\t")
        id = splits[0]
        tuples = [eval(t) for t in splits[1:]]
        results = []
        for tup in tuples:
            (lrmsd, irmsd, fnat) = tup
            results.append((fnat, irmsd, lrmsd))  # Store as (fnat, irmsd, lrmsd)
        res_dict[id] = results
    return res_dict


def is_rigid(key, rigid_models):
    """Checks if model is classified as a hard or rigid model. Returns True if rigid, False otherwise.
    
    Args:
        key (str): Model ID
        rigid_models (list): List of rigid models
    
    return:
        bool: True if rigid, False otherwise
    """
    return key in rigid_models


def main():
    sampling_file = argv[1]  # Path to your sampling results file
    res_dict = parse_sampling_results(sampling_file)

    # path to output file
    outfilename = Path(argv[2])
    # specify the figure size
    plt.figure(figsize=(10, 6))
    ax = plt.gca()

    # Set the font size for the plot
    matplotlib.rcParams.update({'font.size': 18})
    stack_h = 1
    max_stack_v = 0
    stack_h_labels = []
    
    # Define the legend data for the rigid models
    legend_data = [
        (plt.Rectangle((0, 0), 1, 1, fc="darkgreen"), 'High'),
        (plt.Rectangle((0, 0), 1, 1, fc="green"), 'Medium'),
        (plt.Rectangle((0, 0), 1, 1, fc="lightgreen"), 'Acceptable'),
        (plt.Rectangle((0, 0), 1, 1, fc="gray"), 'Incorrect')
    ]
    legend_proxies, legend_labels = zip(*legend_data)
    
    # Define the legend data for medium models
    legend_data_medium = [
        (plt.Rectangle((0, 0), 1, 1, fc="darkblue"), 'High'),
        (plt.Rectangle((0, 0), 1, 1, fc="blue"), 'Medium'),
        (plt.Rectangle((0, 0), 1, 1, fc="lightblue"), 'Acceptable'),
        (plt.Rectangle((0, 0), 1, 1, fc="gray"), 'Incorrect')
    ]
    legend_proxies_medium, legend_labels_medium = zip(*legend_data_medium)

    # Sort the models into rigid and medium models (This is now hardcoded, but I think it is better if it isn't)
    rigid_models = sorted(['1ao7', '1mwa', '2bnr', '2nx5', '2pye', '3dxa', '3pwp','3qdg', '3qdj', '3utt', '3vxr', '3vxs', '5c0a', '5c0b', '5c0c', '5c07', '5c09', '5hyj', '5ivx', '5nme', '5nmf'])
    medium_models = sorted([model for model in res_dict.keys() if model not in rigid_models])
    stack_h_labels = rigid_models + medium_models

    for key in rigid_models:
        rmsds = res_dict[key]
        rigid = True
        stack_v = 0
        for rmsd in rmsds:
            # Check the quality of the model and color code it
            stack_color = _color_code(rmsd, rigid)
            if stack_color == 'white':
                stack_color = 'gray'
            ax.bar(stack_h, 1, bottom=stack_v, color=stack_color,
                  width=1, align='center', edgecolor='none', linewidth=0)
            stack_v += 1
        max_stack_v = max(max_stack_v, stack_v)
        stack_h += 2

    for key in medium_models:
        rmsds = res_dict[key]
        rigid = False
        stack_v = 0
        for rmsd in rmsds:
            # Check the quality of the model and color code it
            stack_color = _color_code(rmsd, rigid)
            if stack_color == 'white':
                stack_color = 'gray'
            ax.bar(stack_h, 1, bottom=stack_v, color=stack_color,
                  width=1, align='center', edgecolor='none', linewidth=0)
            stack_v += 1
        max_stack_v = max(max_stack_v, stack_v)
        stack_h += 2

    # Set the x-axis limits
    ax.set_xlim((0, stack_h - 1))
    ax.set_ylim((-1, max_stack_v + 1))

    # Set the y-axis label and the fontsize
    ax.set_ylabel('Energy rank \n(lower is better)', fontsize=14)  

    # Update for uniform x-axis label spacing with monospaced font and increased size
    stack_h_labels_pos = np.arange(1, stack_h, 2)  # Ensure uniform positioning
    ax.set_xticks(stack_h_labels_pos)
    ax.set_xticklabels(stack_h_labels, rotation=90, fontsize=12, fontfamily='monospace')  # Set monospaced font and increased size

    ax.xaxis.set_ticks_position('none')
    
    # Adjust plot position to accommodate longer x-axis labels
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15,
                     box.width, box.height * 0.82])

    legend_left = ax.legend(legend_proxies, legend_labels,
                            loc='center left', bbox_to_anchor=(0, -0.22), ncol=2,
                            fontsize=8)

    legend_right = ax.legend(legend_proxies_medium, legend_labels_medium,
                             loc='center right', bbox_to_anchor=(1, -0.22), ncol=2,
                             fontsize=8)
    ax.add_artist(legend_left)
    ax.add_artist(legend_right)

    # Set the title and the fontsize
    plt.yticks(fontsize=12)
    plt.title("Melquiplot benchmark", fontsize=18)
    plt.savefig(outfilename, format='pdf', dpi=1500)

if __name__ == "__main__":
    main()