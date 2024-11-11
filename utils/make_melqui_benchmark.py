"""
Name: make_melqui_benchmark.py
Function: This script is used to generate a melquiplot for benchmarking results. This script makes use of hardcoded data for rigid and non-rigid models. The color difference is determined by that data.
Date: 11-11-2024
Author: Nils Smit, Yannick Aarts
"""

"""
Example usage: python make_melqui_benchmark.py data.tsv output
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
import numpy as np
import pandas as pd
from sys import argv
from pathlib import Path

# List of rigid models
rigid_models = [
    '1ao7', '1mwa', '2bnr', '2nx5', '2pye', '3dxa', '3pwp', '3qdg', '3qdj', '3utt', 
    '3vxr', '3vxs', '5c0a', '5c0b', '5c0c', '5c07', '5c09', '5hyj', '5ivx', '5nme', '5nmf'
]

# Define the color mapping based on model type (rigid or non-rigid) with different shades of blue
model_color_map = {
    # easy cases that have minimal conformational changes
    "rigid": {
        "high": "#003600",  # Dark green for high-quality rigid
        "medium": "green",          # Green for medium-quality rigid
        "acceptable": "lightgreen", # Light green for acceptable rigid
        "incorrect": "grey"         # Grey for incorrect rigid models
    },
    # hard cases that have more conformational changes
    "non_rigid": {
        "high": "midnightblue",  # Dark blue for high-quality non-rigid
        "medium": "blue",        # Blue for medium-quality non-rigid
        "acceptable": "lightblue",    # Light blue for acceptable non-rigid
        "incorrect": "grey"    # grey for incorrect non-rigid models
    }
}

def read_tsv_file(tsv_file):
    """Read TSV file and return DataFrame.
    
    Args:
        tsv_file (str): Path to the TSV file
    
    Returns:
        pd.DataFrame: DataFrame containing the data from the TSV file
    """
    return pd.read_csv(tsv_file, sep='\t', index_col=0)

def create_melqui_plot(tsv_file, outfilename):
    """Create a melquiplot based on the input TSV file.
    """
    # Read TSV data
    df = read_tsv_file(tsv_file)

    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 6))  # Adjusted to fit both rigid and non-rigid models
    matplotlib.rcParams.update({'font.size': 12})

    # Determine the number of rigid and non-rigid models
    rigid_df = df.loc[:, df.columns.isin(rigid_models)]
    non_rigid_df = df.loc[:, ~df.columns.isin(rigid_models)]

    # Loop through the rigid models and plot them
    for col_idx, model in enumerate(rigid_df.columns):
        for row_idx, quality in enumerate(rigid_df[model]):
            # Rigid models color based on quality
            color = model_color_map["rigid"].get(quality, "white")
            
            # Plot each rank as a bar with appropriate color for rigid models
            ax.bar(col_idx * 2 + 1, 1, bottom=row_idx, color=color, width=1, 
                   align='center', edgecolor='none', linewidth=0)

    # Loop through the non-rigid models and plot them
    for col_idx, model in enumerate(non_rigid_df.columns):
        for row_idx, quality in enumerate(non_rigid_df[model]):
            # Non-rigid models color based on quality
            color = model_color_map["non_rigid"].get(quality, "white")
            
            # Plot each rank as a bar with appropriate color for non-rigid models
            ax.bar((col_idx + len(rigid_df.columns)) * 2 + 1, 1, bottom=row_idx, color=color, width=1, 
                   align='center', edgecolor='none', linewidth=0)

    # Set the axis limits
    ax.set_xlim((0, (len(rigid_df.columns) + len(non_rigid_df.columns)) * 2))
    ax.set_ylim((0, len(df.index) + 1))
    
    # X-axis label configuration
    ax.set_xticks(np.arange(1, (len(rigid_df.columns) + len(non_rigid_df.columns)) * 2, 2))
    ax.set_xticklabels(list(rigid_df.columns) + list(non_rigid_df.columns), rotation=90, 
                       fontfamily='monospace', fontsize=10)

    # Adjust y-axis and add label
    ax.set_ylabel('Energy rank \n(lower is better)', fontsize=14)

    # Set title with additional space above
    plt.title("Melquiplot", fontsize=18)

    # Create legend for rigid and non-rigid model types with quality color mapping
    legend_patches = [
        mpatches.Patch(color=model_color_map["rigid"]["high"], label="Rigid: High Quality"),
        mpatches.Patch(color=model_color_map["rigid"]["medium"], label="Rigid: Medium"),
        mpatches.Patch(color=model_color_map["rigid"]["acceptable"], label="Rigid: Acceptable"),
        mpatches.Patch(color=model_color_map["rigid"]["incorrect"], label="Rigid: Incorrect"),
        mpatches.Patch(color=model_color_map["non_rigid"]["high"], label="Non-rigid: High Quality"),
        mpatches.Patch(color=model_color_map["non_rigid"]["medium"], label="Non-rigid: Medium"),
        mpatches.Patch(color=model_color_map["non_rigid"]["acceptable"], label="Non-rigid: Acceptable"),
        mpatches.Patch(color=model_color_map["non_rigid"]["incorrect"], label="Non-rigid: Incorrect")
    ]
    
    plt.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.15), 
               title="Model Quality", fontsize=10, title_fontsize='12', ncol=4)

    # Adjust layout to give more space to the y-axis label and the title
    plt.subplots_adjust(bottom=0.25, top=0.9, left=0.15, right=0.9)

    # Save the plot
    plt.savefig(str(outfilename) + '.pdf', format='pdf', dpi=1500)
    print(f"Plot saved to {outfilename}")

if __name__ == "__main__":
    tsv_file = argv[1]
    outfilename = Path(argv[2])
    create_melqui_plot(tsv_file, outfilename)
