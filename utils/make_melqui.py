# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:47:09 2024

@author: yanic
Replicate melqui plot
"""
#imports
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
from pathlib import Path
import numpy as np
import matplotlib.colors as mcolors

def check_fnat(fnat):
    """Function to check if fnat is high, medium, acceptable, incorrect coded as 3, 2, 1, 0

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

#functions
def _color_code(value, rigid):  # a tuple (fnat, irmsd, lrmsd)
    (fnat, irmsd, lrmsd) = value
    fnat_score = check_fnat(fnat)
    irmsd_score = check_irmsd(irmsd)
    lrmsd_score = check_lrmsd(lrmsd)
    if rigid:
        if fnat_score == 3 and (irmsd_score == 3 or lrmsd_score == 3):
            return '#003600'  # 'darkgreen'
        elif fnat_score >= 2 and (irmsd_score >= 2 or lrmsd_score >= 2):
            return 'green'
        elif fnat_score >= 1 and (irmsd_score >= 1 or lrmsd_score >= 1):
            return 'lightgreen'
        else:
            return 'white'  # 'gray'
    else:
        if fnat_score == 3 and (irmsd_score == 3 or lrmsd_score == 3):
            return 'midnightblue'  # 'darkblue'
        elif fnat_score >= 2 and (irmsd_score >= 2 or lrmsd_score >= 2):
            return 'blue'
        elif fnat_score >= 1 and (irmsd_score >= 1 or lrmsd_score >= 1):
            return 'lightblue'
        else:
            return 'white'  # 'gray'

def parse_sampling_results(sampling_file):
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
    """
    return key in rigid_models

def main():
    rmsd_types = ["lrmsd", "irmsd", "fnat"]
    sampling_file = r'results_swifttcr/' + 'combined_all_rmsd_3'  # Path to your sampling results file
    res_dict = parse_sampling_results(sampling_file)

    outfilename = Path(r"results_swifttcr/" + "melquiplot" + ".pdf")
    f = plt.figure(figsize=(10, 6))  # Specify figure size for better space distribution
    ax = plt.gca()

    matplotlib.rcParams.update({'font.size': 18})
    stack_h = 1
    max_stack_v = 0
    stack_h_labels = []
    legend_data = [
        (plt.Rectangle((0, 0), 1, 1, fc="darkgreen"), 'High'),
        (plt.Rectangle((0, 0), 1, 1, fc="green"), 'Medium'),
        (plt.Rectangle((0, 0), 1, 1, fc="lightgreen"), 'Acceptable'),
        (plt.Rectangle((0, 0), 1, 1, fc="gray"), 'Incorrect')
    ]
    legend_proxies, legend_labels = zip(*legend_data)
    legend_data_medium = [
        (plt.Rectangle((0, 0), 1, 1, fc="darkblue"), 'High'),
        (plt.Rectangle((0, 0), 1, 1, fc="blue"), 'Medium'),
        (plt.Rectangle((0, 0), 1, 1, fc="lightblue"), 'Acceptable'),
        (plt.Rectangle((0, 0), 1, 1, fc="gray"), 'Incorrect')
    ]
    legend_proxies_medium, legend_labels_medium = zip(*legend_data_medium)

    rigid_models = sorted(['1ao7', '1mwa', '2bnr', '2nx5', '2pye', '3dxa', '3pwp','3qdg', '3qdj', '3utt', '3vxr', '3vxs', '5c0a', '5c0b', '5c0c', '5c07', '5c09', '5hyj', '5ivx', '5nme', '5nmf'])  # Update this list as needed
    medium_models = sorted([model for model in res_dict.keys() if model not in rigid_models])
    stack_h_labels = rigid_models + medium_models

    for key in rigid_models:
        rmsds = res_dict[key]
        rigid = True
        stack_v = 0
        for rmsd in rmsds:
            stack_color = _color_code(rmsd, rigid)  # Pass the entire rmsd tuple
            if stack_color == 'white':
                stack_color = 'gray'
            series = ax.bar(stack_h, 1, bottom=stack_v, color=stack_color,
                            width=1, align='center', edgecolor='none', linewidth=0)
            stack_v += 1
        max_stack_v = max(max_stack_v, stack_v)
        stack_h += 2  # Separated by an empty series

    for key in medium_models:
        rmsds = res_dict[key]
        rigid = False
        stack_v = 0
        for rmsd in rmsds:
            stack_color = _color_code(rmsd, rigid)  # Pass the entire rmsd tuple
            if stack_color == 'white':
                stack_color = 'gray'
            series = ax.bar(stack_h, 1, bottom=stack_v, color=stack_color,
                            width=1, align='center', edgecolor='none', linewidth=0)
            stack_v += 1
        max_stack_v = max(max_stack_v, stack_v)
        stack_h += 2  # Separated by an empty series

    # Aesthetics
    ax.set_xlim((0, stack_h - 1))  # Centered
    ax.set_ylim((-1, max_stack_v + 1))

    ax.set_ylabel('Energy rank \n(lower is better)', fontsize=12)  # Adjusted font size

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

    plt.yticks(fontsize=8)
    plt.title("Melquiplot benchmark", fontsize=16)  # Adjusted title size for better balance
    plt.savefig(outfilename, format='pdf', dpi=1500)

if __name__ == "__main__":
    main()