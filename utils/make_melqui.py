"""
Name: make_melqui.py
Function: This script is used to generate a melquiplot for results. The format of the input file is a TSV file with models as columns and quality levels as rows. 
Date: 11-11-2024
Author: Nils Smit, Yannick Aarts
"""

"""
Example usage: python make_melqui.py data.tsv output '{"high": "#003600", "medium": "green", "acceptable": "lightgreen", "incorrect": "grey"}'

If the quality_color_map argument is not provided, the default color map will be used.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
import numpy as np
import pandas as pd
from sys import argv
from pathlib import Path
import json

def read_tsv_file(tsv_file):
    """Read TSV file and return DataFrame."""
    return pd.read_csv(tsv_file, sep='\t', index_col=0)


def main():
    # File paths
    tsv_file = argv[1]
    outfilename = Path(argv[2])
    
    # Default color map
    default_quality_color_map = {
        "high": "#003600",
        "medium": "green",
        "acceptable": "lightgreen",
        "incorrect": "grey"
    }
    
    # Parse quality_color_map argument or use default
    try:
        quality_color_map = json.loads(argv[3])
    except (IndexError, json.JSONDecodeError):
        quality_color_map = default_quality_color_map

    # Read TSV data
    df = read_tsv_file(tsv_file)

    # Set up the plot
    fig, ax = plt.subplots(figsize=(10, 6))  # Create a figure and axis
    matplotlib.rcParams.update({'font.size': 12})
    
    # Loop through the DataFrame and plot each cell with the appropriate color
    for col_idx, model in enumerate(df.columns):
        for row_idx, quality in enumerate(df[model]):
            # Determine the color based on the quality value
            color = quality_color_map.get(quality, "white")  # Default to white if quality is not recognized

            # Plot each rank as a bar with appropriate color
            ax.bar(col_idx * 2 + 1, 1, bottom=row_idx, color=color, width=1, 
                   align='center', edgecolor='none', linewidth=0)

    # Set the axis limits
    ax.set_xlim((0, len(df.columns) * 2))
    ax.set_ylim((0, len(df.index) + 1))
    
    # X-axis label configuration
    ax.set_xticks(np.arange(1, len(df.columns) * 2, 2))
    ax.set_xticklabels(df.columns, rotation=90, fontfamily='monospace', fontsize=10)

    # Adjust y-axis and add label
    ax.set_ylabel('Energy rank \n(lower is better)', fontsize=14)

    # Set title with additional space above
    plt.title("Melquiplot", fontsize=18)

    # Create legend from quality_color_map and position it below the plot
    legend_patches = [mpatches.Patch(color=color, label=quality) 
                      for quality, color in quality_color_map.items()]
    plt.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.15), 
               title="Quality Levels", fontsize=10, title_fontsize='12', ncol=4)

    # Adjust layout to give more space to the y-axis label and the title
    plt.subplots_adjust(bottom=0.25, top=0.9, left=0.15, right=0.9)

    # Save the plot
    plt.savefig(str(outfilename) + '.pdf', format='pdf', dpi=1500)
    print(f"Plot saved to {outfilename}")

if __name__ == "__main__":
    main()