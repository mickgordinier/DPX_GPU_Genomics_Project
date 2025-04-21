import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def compare_versions(df, components, title, x_label, y_label, outfile):
    df.fillna(0, inplace=True)
    plt.style.use('ggplot')

    # Bar positions and width
    x = range(len(df))
    bar_width = 0.6

    # Plotting each stack
    bottom = [0] * len(df)

    plt.figure(figsize=(15, 6))

    for comp in components:
        plt.bar(df[x_label], df[comp], bottom=bottom, label=comp)
        bottom = [i + j for i, j in zip(bottom, df[comp])]

    # Annotate bars with average execution time for each kernel
    for idx, version in enumerate(df[x_label]):
        total = sum(df.loc[idx, comp] for comp in components)
        plt.text(idx, total + 0.02, f"{total:.2f} us", ha='center', va='bottom', fontsize=9)

    # Labels and legend
    plt.title(title, fontsize=14)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend(title="Component", loc="upper right")
    plt.tight_layout()
    plt.grid(axis='y', linestyle='--', alpha=0.6)

    plt.savefig(outfile, dpi=300)


version_components =  ['Memory Management', 'Kernel Execution', 'Backtracking', 'Printing', 'Misc']
version_df = pd.read_excel('Timing_Analysis.xlsx', sheet_name='Results', header=1)
compare_versions(version_df, version_components, "Average Execution Time Breakdown per Kernel Version", "Version", "Time (us)", "kernel_timing_breakdown.png")

gpu_comparison_components = ['Memory Management', 'Kernel Execution', 'Printing']
gpu_df = pd.read_excel('Timing_Analysis.xlsx', sheet_name='V12 Comparison', header=1)
compare_versions(gpu_df, gpu_comparison_components, "V12 Execution Time Breakdown Across GPUs", "GPU", "Time (ms)", "v12_gpu_comparison.png")

    