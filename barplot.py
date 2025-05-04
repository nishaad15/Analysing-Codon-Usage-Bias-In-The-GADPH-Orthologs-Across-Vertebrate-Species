import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("codon_usage_table.tsv", sep="\t")

codons = df['Codon']
species = df.columns[1:]
num_species = len(species)


fig, axes = plt.subplots(3, 2, figsize=(16, 12), sharey=True)
axes = axes.flatten() 

for i, sp in enumerate(species):
    ax = axes[i]
    ax.bar(codons, df[sp], color='skyblue', edgecolor='black')
    ax.set_title(sp, fontsize=12)
    ax.set_xticks(np.arange(len(codons)))
    ax.set_xticklabels(codons, rotation=90, fontsize=8)
    if i % 2 == 0: 
        ax.set_ylabel("Codon Frequency")
    ax.set_xlabel("Codon", fontsize=9)

if num_species < 6:
    for j in range(num_species, 6):
        fig.delaxes(axes[j])


plt.suptitle("Codon Usage Comparison Across Species", fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("codon_usage_grid_3x2.png", dpi=300)
plt.show()
