import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_excel('Timing_Analysis.xlsx', sheet_name='Results', header=1)
df.fillna(0, inplace=True)
plt.style.use('ggplot')

# Bar positions and width
x = range(len(df))
bar_width = 0.6

# Plotting each stack
bottom = [0] * len(df)
components = ['Memory Management', 'Kernel Execution', 'Backtracking', 'Printing', 'Misc']

plt.figure(figsize=(15, 6))

for comp in components:
    plt.bar(df['Version'], df[comp], bottom=bottom, label=comp)
    bottom = [i + j for i, j in zip(bottom, df[comp])]

# Annotate bars with average execution time for each kernel
for idx, version in enumerate(df['Version']):
    total = sum(df.loc[idx, comp] for comp in components)
    plt.text(idx, total + 0.02, f"{total:.2f} us", ha='center', va='bottom', fontsize=9)

# Labels and legend
plt.title("Average Execution Time Breakdown per Kernel Version", fontsize=14)
plt.xlabel("Kernel Version")
plt.ylabel("Time (us)")
plt.legend(title="Component", loc="upper right")
plt.tight_layout()
plt.grid(axis='y', linestyle='--', alpha=0.6)

plt.savefig("kernel_timing_breakdown.png", dpi=300)
