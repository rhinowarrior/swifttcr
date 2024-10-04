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
def _color_code(value, rigid):#a tuple (fnat, irmsd, lrmsd)
    (fnat, irmsd, lrmsd) = value
    fnat_score = check_fnat(fnat)
    irmsd_score = check_irmsd(irmsd)
    lrmsd_score = check_lrmsd(lrmsd)
    #print (fnat_score, irmsd_score, lrmsd_score)
    #print(fnat_score == 2 and (irmsd_score == 2 or lrmsd_score == 2))
    if rigid:
        if fnat_score == 3 and (irmsd_score == 3 or lrmsd_score == 3):
            return '#003600'#'darkgreen'
        elif fnat_score >= 2 and (irmsd_score >= 2 or lrmsd_score >= 2):
            return 'green'
        elif fnat_score >= 1 and (irmsd_score >= 1 or lrmsd_score >= 1):
            return 'lightgreen'
        else:
            return 'white'#'gray'
    else:
        if fnat_score == 3 and (irmsd_score == 3 or lrmsd_score == 3):
            return 'midnightblue'#'darkblue'
        elif fnat_score >= 2 and (irmsd_score >= 2 or lrmsd_score >= 2):
            return 'blue'
        elif fnat_score >= 1 and (irmsd_score >= 1 or lrmsd_score >= 1):
            return 'lightblue'
        else:
            return 'white'#'gray'

def read_data(filename):
    fh = open(filename)
    lines = fh.readlines()
    fh.close()
    res_dict = {}
    for line in lines:
        splits = line.split()
        case_id = splits[0]
        rmsd_values = splits[1:]
        rmsd_float = [float(x) for x in rmsd_values]
        data = rmsd_float
        res_dict[case_id] = data
    return res_dict

def merge_type_dict(type_dict):
    merged_dict = {}
    fnat_dict = type_dict["fnat"]
    irmsd_dict = type_dict["irmsd"]
    lrmsd_dict = type_dict["lrmsd"]
    common_keys = fnat_dict.keys() & irmsd_dict.keys() & lrmsd_dict.keys()
    for k in common_keys:
        merged_dict[k] = (zip (fnat_dict[k], irmsd_dict[k], lrmsd_dict[k]))
            
    return merged_dict

def is_rigid(key, rigid_models):
    """Checks if model is clasified as a hard or rigid model. Returns True if rigid, False otherwise.
    """
    if key in rigid_models:
        return True
    else:
        return False

def luminance(color):
    """Calculate the luminance of a color.

    Args:
        color (tuple): A tuple (R, G, B) where each value is between 0 and 255.

    Returns:
        float: The luminance of the color.
    """
    r, g, b = color
    return 0.299 * r + 0.587 * g + 0.114 * b

def get_rgb(color):
    """Convert a color name to an RGB tuple, or return the tuple if already an RGB tuple.

    Args:
        color (str or tuple): The color name or RGB tuple.

    Returns:
        tuple: The RGB tuple corresponding to the color.
    """
    if isinstance(color, str):
        # Convert color name to RGB tuple
        return tuple(int(c * 255) for c in mcolors.to_rgb(color))
    elif isinstance(color, tuple) and len(color) == 3:
        return color
    else:
        raise ValueError("Color must be a color name or an RGB tuple")

def color_max(color1, color2):
    """Return the darkest of the two colors based on luminance.

    Args:
        color1 (str or tuple): The first color (name or RGB tuple).
        color2 (str or tuple): The second color (name or RGB tuple).

    Returns:
        str or tuple: The color with the lower luminance (i.e., the darker color).
    """
    rgb1 = get_rgb(color1)
    rgb2 = get_rgb(color2)
    return color1 if luminance(rgb1) < luminance(rgb2) else color2

def darkest_color(colors):
    """Return the darkest color from a list of colors.

    Args:
        colors (list): A list of colors (names or RGB tuples).

    Returns:
        str or tuple: The darkest color based on luminance.
    """
    if not colors:
        raise ValueError("The color list should not be empty.")
    
    darkest = colors[0]
    darkest_luminance = luminance(get_rgb(darkest))

    for color in colors[1:]:
        current_luminance = luminance(get_rgb(color))
        if current_luminance < darkest_luminance:
            darkest = color
            darkest_luminance = current_luminance
    
    return darkest

def sort_colors_by_darkest(colors):
    """Sort a list of colors with the darkest color first based on luminance.

    Args:
        colors (list): A list of colors (names or RGB tuples).

    Returns:
        list: The sorted list of colors with the darkest color first.
    """
    return sorted(colors, key=lambda color: luminance(get_rgb(color)))


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
            results.append((fnat, irmsd, lrmsd))
        res_dict[id] = results
    return res_dict

        





    


def make_rmsd_plot():
    
    rmsd_types = ["fnat", "irmsd", "lrmsd"]
    type_dict = {}
    for threshold in [2,3,4,5,6,7,8,9]:
        for rmsd_type in rmsd_types:
            filename = r'C:\Users\aarts026\Documents\python_programs\data\jandata\jandata\{}_benchmark_60_{}'.format(rmsd_type, threshold)
            res_type_dict = read_data(filename)
            type_dict[rmsd_type] = res_type_dict
        #outfilename = Path(r"C:\Users\aarts026\Documents\python_programs\images" , "melquiplot_" + str(str(filename).split("_benchmark_")[1]) + ".pdf")
        res_dict = merge_type_dict(type_dict)
        rigid_models = sorted(['1ao7', '1mwa', '2bnr', '2nx5', '2pye', '3dxa', '3pwp', '3qdj', '3qdg', '3utt', '3vxr', '3vxs', '5c0a', '5c0b', '5c0c', '5c07', '5c09', '5hyj', '5ivx', '5nme', '5nmf'])
        medium_models = sorted([model for model in res_dict.keys() if model not in rigid_models])
        labels = [1,5,10,20,50,100]
        # Set square size and spacing
        square_size = 0.8
        spacing = 0.2
        # Create plot
        fig, ax = plt.subplots()
        i = 0
        for model in rigid_models + medium_models:
            rigid = is_rigid(model, rigid_models)
            rmsds = res_dict[model]
            
            j=0
            highest_color = 'white'
            colors = []
            for rmsd in rmsds:
                rect_color = _color_code(rmsd, rigid)
                colors.append(rect_color)
            for label in labels:
                highest_color = darkest_color(colors[:label])
                rect = plt.Rectangle((i - square_size / 2, j - square_size / 2), square_size, square_size, color=highest_color)
                ax.add_patch(rect)
                j+=1
            i+=1

        # Set plot margins and aspect ratio
        plt.margins(0.2)
        plt.gca().set_aspect('equal', adjustable='box')
        
        # Set x and y axis limits
        plt.xlim(-0.5, len(res_dict) - 0.5)
        plt.ylim(-0.5, len(labels) - 0.5)
        
        # Set x and y axis labels
        ax.set_xticks(np.arange(len(res_dict)))
        ax.set_xticklabels(list(map(lambda s: s.upper(), rigid_models + medium_models)), rotation=90, fontsize=8, fontname='DejaVu Sans Mono')
        ax.set_yticks(np.arange(len(labels)) - 0.2)  # Adjust the value as needed
        ax.set_yticklabels(['Top ' + str(label) for label in labels], fontsize=8, fontname='Helvetica' )
        ax.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False)
        ax.set_title("Performance on Top N models",  fontsize=8, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        legend_data = [ (plt.Rectangle((0, 0), 1, 1, fc="darkgreen"), 'High'),
                    (plt.Rectangle((0, 0), 1, 1, fc="green"), 'Medium'),
                    (plt.Rectangle((0, 0), 1, 1, fc="lightgreen"), 'Acceptable'),
                    (plt.Rectangle((0, 0), 1, 1, fc="gray"), 'Incorrect') ]
        legend_proxies, legend_labels = zip(*legend_data)
        legend_data_medium = [ (plt.Rectangle((0, 0), 1, 1, fc="darkblue"), 'High'),
                    (plt.Rectangle((0, 0), 1, 1, fc="blue"), 'Medium'),
                    (plt.Rectangle((0, 0), 1, 1, fc="lightblue"), 'Acceptable'),
                    (plt.Rectangle((0, 0), 1, 1, fc="gray"), 'Incorrect') ]
        legend_proxies_medium, legend_labels_medium = zip(*legend_data_medium)
        legend_left = ax.legend(legend_proxies, legend_labels,
                loc='center left', bbox_to_anchor=(0, -0.22), ncol=2,
                fontsize=8)
        
        legend_right = ax.legend(legend_proxies_medium, legend_labels_medium,
                loc='center right', bbox_to_anchor=(1, -0.22), ncol=2,
                fontsize=8)
        ax.add_artist(legend_left)
        ax.add_artist(legend_right)
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
        outfilename = Path(r"C:\Users\aarts026\Documents\python_programs\images" , "rmsdheatmap_" + str(str(filename).split("_benchmark_")[1]) + ".pdf")
        plt.savefig(outfilename, format='pdf', dpi=1500)



#main
def main():
    rmsd_types = ["fnat", "irmsd", "lrmsd"]
    type_dict = {}
    for threshold in [2,3,4,5,6,7,8,9]:
        #res_dict = parse_sampling_results(r'C:\Users\aarts026\Documents\python_programs\data\results_all.txt')
        for rmsd_type in rmsd_types:
            filename = r'C:\Users\aarts026\Documents\python_programs\data\jandata\jandata\{}_benchmark_60_{}'.format(rmsd_type, threshold)
            res_type_dict = read_data(filename)
            type_dict[rmsd_type] = res_type_dict
        outfilename = Path(r"C:\Users\aarts026\Documents\python_programs\images" , "melquiplot_" + str(str(filename).split("_benchmark_")[1]) + ".pdf")
        res_dict = merge_type_dict(type_dict)
        #outfilename = Path(r"C:\Users\aarts026\Documents\python_programs\images" , "melquiplot_sampling.pdf")
        
        f = plt.figure()#figsize=(270/25.4, 320/25.4), dpi=300)
        ax = plt.gca()
        
        matplotlib.rcParams.update({'font.size': 18})
        stack_h = 1
        max_stack_v = 0
        stack_h_labels = []
        legend_data = [ (plt.Rectangle((0, 0), 1, 1, fc="darkgreen"), 'High'),
                    (plt.Rectangle((0, 0), 1, 1, fc="green"), 'Medium'),
                    (plt.Rectangle((0, 0), 1, 1, fc="lightgreen"), 'Acceptable'),
                    (plt.Rectangle((0, 0), 1, 1, fc="gray"), 'Incorrect') ]
        legend_proxies, legend_labels = zip(*legend_data)
        legend_data_medium = [ (plt.Rectangle((0, 0), 1, 1, fc="darkblue"), 'High'),
                    (plt.Rectangle((0, 0), 1, 1, fc="blue"), 'Medium'),
                    (plt.Rectangle((0, 0), 1, 1, fc="lightblue"), 'Acceptable'),
                    (plt.Rectangle((0, 0), 1, 1, fc="gray"), 'Incorrect') ]
        legend_proxies_medium, legend_labels_medium = zip(*legend_data_medium)
        
        rigid_models = sorted(['1ao7', '1mwa', '2bnr', '2nx5', '2pye', '3dxa', '3pwp', '3qdj', '3utt', '3vxr', '3vxs', '5c0a', '5c0b', '5c0c', '5c07', '5c09', '5hyj', '5ivx', '5nme', '5nmf'])
        medium_models = sorted([model for model in res_dict.keys() if model not in rigid_models])
        stack_h_labels = rigid_models + medium_models
        for key in rigid_models:
            rmsds = res_dict[key]
            rigid = True
            stack_v = 0
            for rmsd in rmsds:
                stack_color = _color_code(rmsd, rigid)#rewrite to take tuple (irmsd, lrmsd, fnat)
                if stack_color == 'white':
                    stack_color = 'gray'
                #series = ax.bar(stack_h, 1, bottom=stack_v, color=stack_color,
                #                width=0.5, align='center', edgecolor='none', linewidth=0)
                series = ax.bar(stack_h, 1, bottom=stack_v, color=stack_color,
                                width=1, align='center', edgecolor='none', linewidth=0)
                stack_v += 1
            max_stack_v = max(max_stack_v, stack_v)
            stack_h += 2 # Separated by an empty series
        for key in medium_models:
            rmsds = res_dict[key]
            rigid = False
            stack_v = 0
            for rmsd in rmsds:
                stack_color = _color_code(rmsd, rigid)#rewrite to take tuple (irmsd, lrmsd, fnat)
                if stack_color == 'white':
                    stack_color = 'gray'
                #series = ax.bar(stack_h, 1, bottom=stack_v, color=stack_color,
                #                width=0.5, align='center', edgecolor='none', linewidth=0)
                series = ax.bar(stack_h, 1, bottom=stack_v, color=stack_color,
                                width=1, align='center', edgecolor='none', linewidth=0)
                stack_v += 1
            max_stack_v = max(max_stack_v, stack_v)
            stack_h += 2 # Separated by an empty series
        
        # Aesthetics
        ax.set_xlim((0, stack_h-1)) # Centered
        ax.set_ylim((-1, max_stack_v + 1))

        ax.set_ylabel('Energy rank \n(lower is better)', fontsize=8)
        #ax.set_title(os.path.basename(args.data_file))

        stack_h_labels_pos = range(1, stack_h, 2)
        ax.set_xticks(stack_h_labels_pos)
        ax.set_xticklabels(stack_h_labels, rotation=90, fontsize = 8)
        ax.xaxis.set_ticks_position('none')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.12,
                        box.width, box.height * 0.82])

        # Put a legend below current axis
    #    ax.legend(legend_proxies, legend_labels,
    #              loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2,
    #              fontsize='small')
        legend_left = ax.legend(legend_proxies, legend_labels,
                loc='center left', bbox_to_anchor=(0, -0.22), ncol=2,
                fontsize=8)
        
        legend_right = ax.legend(legend_proxies_medium, legend_labels_medium,
                loc='center right', bbox_to_anchor=(1, -0.22), ncol=2,
                fontsize=8)
        ax.add_artist(legend_left)
        ax.add_artist(legend_right)
        # Save figure
        #plt.savefig(args.output, format='eps')
        plt.yticks(fontsize=8)
        plt.title("Melquiplot cluster irmsd < {} Å".format(threshold))
        plt.savefig(outfilename, format='pdf', dpi=1500)
    #make_rmsd_plot()


if __name__ == "__main__":
    main()
