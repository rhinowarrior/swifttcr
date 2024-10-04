# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:52:30 2023

@author: yanic
Script to plot T1, T5, T10, T20, T50, T100 for every case, colored by Capri criteria.

"""
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib import ticker
from matplotlib.patches import Circle, Rectangle
from matplotlib import rc

import matplotlib.pyplot as plt
import numpy as np
import json



import pandas as pd

import plotly.express as px
import plotly.graph_objs as go

n_samples = ''
lrmsd_f = r'/home/nils/swifttcr/lrmsd_rmsd3'
irmsd_f = r'/home/nils/swifttcr/irmsd_rmsd3'
fnat_f  = r'/home/nils/swifttcr/fnat_rmsd3'
succes_plot_f = 'PIPER_success_plot_'+n_samples
success_plot_f = 'PIPER_success_rate_plot_'+n_samples


def main():

    lrmsd_dict = parse_results(lrmsd_f)
    irmsd_dict = parse_results(irmsd_f)
    fnat_dict = parse_results(fnat_f)
    rmsd_plot_name = succes_plot_f
    success_plot_name = success_plot_f

    labels = [1,5,10,20,50,100]

    color_matrix_all = eval_model_qualities(lrmsd_dict, labels,  irmsd_dict, fnat_dict, rmsd_plot_name)
   

    #print(color_matrix_all)
    with open('data.txt', 'w') as f:
        for model in color_matrix_all:
            f.write(model + "\n")
            f.write(str(color_matrix_all[model]) + "\n")

    # # Save to a text file
    # with open('data.txt', 'w') as f:
    #     for model in color_matrix_all:
    #         print(model)
    #     json.dump(color_matrix_all, f)

    # # Convert numpy arrays to lists
    # data = {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in color_matrix_all.items()}

    # # Save to a text file
    # with open('data.txt', 'w') as f:
    #     json.dump(data, f)

    color_matrix_top = calc_max_values([1, 5, 10, 20, 50, 100], color_matrix_all)

    #print(color_matrix_top)
   
    color_matrix = make_rmsd_plot_all_at_one(color_matrix_top, labels, rmsd_plot_name)

    color_matrix[color_matrix == 'blue'] = 'green' 
    color_matrix[color_matrix == 'lightblue'] = 'lightgreen' 
    color_matrix[color_matrix == 'darkblue'] = 'darkgreen' 

    colors = ['darkgreen', 'green', 'lightgreen','white']

    color_counts, percenteages = calculate_color_percentages(color_matrix, labels)
    #print(percenteages)
    plot_success_rate(percenteages, success_plot_name)



def parse_results(datapath):
    """
    Parses a results file and returns a dictionary containing the data.

    Args:
        datapath: The path to the results file.

    Returns:
        A dictionary where the keys are model names and the values are lists of lrmsd values.
    """
    data_dict = {}
    f = open(datapath)
    lines = f.readlines()
    for line in lines:
        splits = line.split()
        model_name = splits[0].strip()
        lrmsd_values = [float(i) for i in splits[1:]]
        data_dict[model_name] = lrmsd_values
    return data_dict

def calc_max_values(indices, data):
    """
    Calculates the maximum values for specified indices from a given data dictionary.

    Args:
        indices (list): A list of indices specifying the number of maximum values to consider.
        data (dict): A dictionary containing the data.

    Returns:
        dict: A dictionary where the keys are model names and the values are lists of maximum values.
    """
    result = {}
    for key, values in data.items():
        if key not in result:
            result[key] = []
        if len(values)>0:
            for index in indices:
                result[key].append(max(values[:index]))
    return result
    
def eval_model_qualities(l_dict, labels ,  i_dict=False, f_dict=False, rmsd_plot_name="rmsd_plot"):

    l_thresholds = [0, 1, 5, 10]  # Threshold values for color levels
    i_thresholds = [0, 1, 2, 4]  # Threshold values for color levels
    f_thresholds = [1, 0.5, 0.3, 0.1]  # Threshold values for color levels
    quality = [3,2,1]

    color_matrix_all ={}
    model_quality = np.empty((37, 100),dtype=int)
    i=0
    for pdb_id, lrmsd_values in l_dict.items():
        
        irmsd_values = i_dict[pdb_id]
        frmsd_values = f_dict[pdb_id]

        print(lrmsd_values,irmsd_values,frmsd_values)

        model_quality = np.zeros((len(f_dict[pdb_id])),dtype=int)
        #print(pdb_id, model_quality)
        for j, (lrmsd, irmsd, frmsd) in enumerate(zip(lrmsd_values, irmsd_values, frmsd_values)):
            
            # for k in range(len(l_thresholds)-1):
                if  ((frmsd >= 0.5) and ((lrmsd <= 1) or  ( irmsd <= 1))):
                    model_quality[j] = 3
                    
                elif  ((frmsd >= 0.3) and ((lrmsd <= 5) or  ( irmsd <= 2))):
                    model_quality[j] = 2
                    
                elif  ((frmsd >= 0.1) and ((lrmsd <= 10) or  (irmsd <= 4))):
                    model_quality[j] = 1
                    
                elif  (frmsd < 0.1 ) :
                    model_quality[j] = 0
                    
                else:
                    print(pdb_id,lrmsd, irmsd, frmsd)

                #print(pdb_id, model_quality)
                color_matrix_all[pdb_id]= model_quality   
                i+=1
    
    return color_matrix_all


def make_rmsd_plot_all_at_one(l_dict, labels, rmsd_plot_name="rmsd_plot"):

    rigid_models = ['1ao7', '1mwa', '2bnr', '2nx5', '2pye', '3dxa', '3pwp','3qdg', '3qdj', '3utt', '3vxr', '3vxs', '5c0a', '5c0b', '5c0c', '5c07', '5c09', '5hyj', '5ivx', '5nme', '5nmf']
    medium_models = [model for model in l_dict if model not in rigid_models]
    # Define color thresholds
    # rigid_colors = ['darkgreen', 'green', 'lightgreen']
    # medium_colors = ['darkblue', 'blue', 'lightblue']

        
    rigid_colors = ['white','lightgreen', 'green', 'darkgreen']
    medium_colors = ['white','lightblue', 'blue', 'darkblue']
    l_thresholds = [0, 1, 5, 10]  # Threshold values for color levels
    i_thresholds = [0, 1, 2, 4]  # Threshold values for color levels
    f_thresholds = [1, 0.5, 0.3, 0.1]  # Threshold values for color levels
    
    # Prepare x and y data
    x = np.arange(len(l_dict))
    y = np.arange(len(labels))
    
    # Set square size and spacing
    square_size = 0.8
    spacing = 0.2
    
    # Create plot
    fig, ax = plt.subplots()

    # Plot rigid models
    model_num = len(rigid_models+medium_models)
    colors = np.zeros((model_num, 6),dtype=object)
    # color_dict = {}
    for i, model in enumerate(sorted(rigid_models) + sorted(medium_models)):

        rmsd_values = l_dict[model]
        print(model, rmsd_values)
        for j,rmsd in enumerate(rmsd_values):
                color=None
                if model in rigid_models:
                            color = rigid_colors[rmsd]
                            
                if model in medium_models:
                            color = medium_colors[rmsd]


                if color:
                    colors[i,j] = color
                    if color =='darkgreen':
                        # color = '#014501'
                        color = '#003600'
                        
                    if color =='darkblue':
                        color = 'midnightblue'
                    rect = plt.Rectangle((i - square_size / 2, j - square_size / 2), square_size, square_size, color=color)
                    ax.add_patch(rect)  
                    

    #   color_dict[model] = colors[i,:]
    # rigid_models[1]='*1mwa'
    # rigid_models[-3]='*5ivx'

    colors[colors == 0] = 'white'

    rigid_models_u = [pdb_id.upper() for pdb_id in rigid_models]
    medium_models_u = [pdb_id.upper() for pdb_id in medium_models]



    # Set plot margins and aspect ratio
    plt.margins(0.2)
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Set x and y axis limits
    plt.xlim(-0.5, len(l_dict) - 0.5)
    plt.ylim(-0.5, len(labels) - 0.5)
    
    # Set x and y axis labels
    ax.set_xticks(np.arange(len(l_dict)))
    #ax.set_xticklabels(t_dict.keys())
    #ax.set_xticklabels(sorted(rigid_models) + sorted(medium_models), rotation=90, fontsize=22, fontname='Arial')  # Rotate and sort labels
    #ax.set_yticklabels(['Top ' + str(label) for label in labels], fontsize=24, fontname='Arial', fontweight='bold')

    ax.set_xticklabels(sorted(rigid_models_u) + sorted(medium_models_u), rotation=90, fontsize=8, fontname='DejaVu Sans Mono')  # Rotate and sort labels
    #ax.set_yticks(np.arange(len(labels)))
    #ax.set_yticklabels(['Top ' + str(label) for label in labels], fontsize=24, fontname='Helvetica', ha='left')
    # Set y-axis tick labels with left alignment

    # Set y-axis tick labels with left alignment
    ax.set_yticks(np.arange(len(labels)) - 0.2)  # Adjust the value as needed
    ax.set_yticklabels(['Top ' + str(label) for label in labels], fontsize=8, fontname='Helvetica' )
    ax.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False)

    # Manually adjust the y-axis tick positions
   
    # ax.set_yticklabels(['Top ' + str(label) for label in labels], fontsize=24, fontname='Helvetica')
    # ax.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, pad=-20)  # Adjust pad as needed



    #ax.yaxis.set_label_coords(-5, 0.5)  # Adjust the coordinates as needed

    #ax.set_yticklabels(['Top ' + str(label) for label in labels],fontsize = 28, fontname='Helvetica')
    
    # Set plot title and labels
    ax.set_title("Performance on Top N models",  fontsize=8, fontweight='bold')
    #ax.set_xlabel("Models")
    #ax.set_ylabel("T")
    
    # Remove ticks and labels on the top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add whitespace between data and x and y axes
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
    
    # Show the plot
    #plt.show()
    plt.savefig(rmsd_plot_name+'.png',bbox_inches='tight',  dpi=300)
    plt.savefig(rmsd_plot_name+'.pdf',bbox_inches='tight',  dpi=300)

    return colors


def calculate_color_percentages(color_matrix, labels):
    # num = np.zeros(len( ))
    # print(num)
    percentages = {}
    color_counts = np.zeros((6,4),dtype=float)
    for i, label in enumerate(labels):

        counts = {'darkgreen': 0, 'green': 0, 'lightgreen': 0, 'white': 0 }

        # counts = {'white': 0, 'lightgreen': 0, 'green': 0, 'darkgreen': 0}
        for model_colors in color_matrix:
            counts[model_colors[i]] += 1
        total = sum(counts.values())
        percentages[label] = {color: (count / total) * 100 for color, count in counts.items()}
        #print(percentages.values())
        color_counts[i]=list(percentages[label].values())
    return color_counts, percentages



def plot_success_rate(data, success_plot_file):
    # Function to reorder dictionary keys
    def reorder_dict_keys(dictionary, keys_order):
        return {key: dictionary[key] for key in keys_order}

    # Reorder the dictionary keys for each row
    for key in data:
        data[key] = reorder_dict_keys(data[key], ['darkgreen', 'green', 'lightgreen', 'white'])

    # Convert the dictionary to DataFrame
    df = pd.DataFrame(data).T.reset_index().rename(columns={'index': 'Top'})
    
    # Melt the DataFrame
    df = df.melt(id_vars=['Top'], var_name='color', value_name='percentage')

    # Convert 'Top' column values to strings
    df['Top'] = df['Top'].astype(str)

    # Reorder the DataFrame columns
    # desired_order = ['darkgreen', 'green', 'lightgreen', 'white']
    desired_order = ['#003600', 'green', 'lightgreen', 'white']
    df['color'] = df['color'].replace('darkgreen', '#003600')
    df['color'] = pd.Categorical(df['color'], categories=desired_order, ordered=True)


    # Create plot with color mapping
    fig = px.bar(df, x="Top", y="percentage", color="color", title="Success Rate",
                 color_discrete_map={color: color for color in desired_order})

    # Update layout to adjust bar spacing
    # fig.update_layout(bargap=0.1, font=dict(size=40),xaxis_title="Top")
    # fig.update_layout(xaxis=dict(tickmode='array', tickvals=list(range(len(df))), ticktext=['Top '+str(val) for val in df['Top']]))

    #fig.update_xaxes(showticklabels=False)
    fig.update_layout(bargap=0.1, font=dict(size=60),xaxis=dict(tickmode='array', ticktext=['Top '+str(val) for val in df['Top']]), title_x=0.5, title=dict(y=0.98), legend_title_text='Model Quality')
    #fig.update_xaxes(showticklabels=True)

    # Update legend labels
    legend_labels = {'#003600': 'High', 'green': 'Medium', 'lightgreen': 'Acceptable', 'white': 'Incorrect'}

    # for color, label in legend_labels.items():
    #     fig.for_each_trace(lambda t: t.update(marker=dict(color=color, size=20)))

    # fig.for_each_trace(lambda t: t.update(name = legend_labels[t.name],
    #                                   legendgroup = legend_labels[t.name],
    #                                   hovertemplate = t.hovertemplate.replace(t.name, legend_labels[t.name])
    #                                  )
    #               )
    # for color, label in legend_labels.items():
    #     fig.for_each_trace(lambda t: t.update(name=label,
    #                                         legendgroup=label,
    #                                         hovertemplate=t.hovertemplate.replace(t.name, label),
    #                                         base=0)
    #                     )
    
    # Update legend labels for each trace
    for color, label in legend_labels.items():
        for trace in fig.data:
            if trace.marker.color == color:
                trace.name = label

    # Update legend group and hover template
    for trace in fig.data:
        label = trace.name
        trace.legendgroup = label
        trace.hovertemplate = trace.hovertemplate.replace(trace.name, label)

    # Add hover text for each trace
    for trace in fig.data:
        trace.hovertext = [f"{label}: {y:.2f}%" for label, y in zip(df['Top'], trace.y)]

    # # Add stacked labels
    # for trace in fig.data:
    #     y_values = trace.y
    #     text_values = [f"{y:.2f}%" for label, y in zip(df['Top'], y_values)]
    #     trace.text = text_values
    #     trace.textposition = 'inside'
    #     trace.textangle = 0  # Set text angle to 0 for horizontal labels



    # Save the plot as an HTML file
    fig.write_html(success_plot_file+".html")



if __name__ == "__main__":
    main()
