import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def compare_versions(df, components, title, x_label, y_label, outfile):
    df.fillna(0, inplace=True)
    plt.style.use('ggplot')

    # Convert microseconds to milliseconds for plotting
    df_ms = df.copy()
    for comp in components:
        df_ms[comp] = df[comp] / 1000

    # Bar positions and width
    x = range(len(df_ms))
    bar_width = 0.6

    # Plotting each stack
    bottom = [0] * len(df_ms)

    plt.figure(figsize=(15, 6))

    for comp in components:
        plt.bar(df_ms[x_label], df_ms[comp], bottom=bottom, label=comp)
        bottom = [i + j for i, j in zip(bottom, df_ms[comp])]

    # Annotate bars with average execution time for each kernel
    for idx, version in enumerate(df_ms[x_label]):
        total = sum(df_ms.loc[idx, comp] for comp in components)
        plt.text(idx, total + 0.02, f"{int(round(total))} ms", 
                ha='center', va='bottom', fontsize=12)

    # Labels and legend - all set to fontsize 12
    plt.title(title, fontsize=12)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    
    # Get current axis and set tick label size
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Legend with fontsize 12
    plt.legend(title="Component", loc="upper right", fontsize=12, title_fontsize=12)
    
    plt.tight_layout()
    plt.grid(axis='y', linestyle='--', alpha=0.6)

    plt.savefig(outfile, dpi=300)


version_components =  ['Memory Management', 'Kernel Execution', 'Backtracking', 'Printing', 'Misc']
version_df = pd.read_excel('Timing_Analysis.xlsx', sheet_name='Results', header=1)
compare_versions(version_df, version_components, "Average Execution Time Breakdown per Kernel Version", "Version", "Time (ms)", "kernel_timing_breakdown.png")

gpu_comparison_components = ['Memory Management', 'Kernel Execution', 'Other']
gpu_df = pd.read_excel('Timing_Analysis.xlsx', sheet_name='V12 Comparison', header=1)
compare_versions(gpu_df, gpu_comparison_components, "V12 Execution Time Breakdown Across GPUs", "GPU", "Time (ms)", "v12_gpu_comparison.png")