import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv("rscu_table_wide.tsv", sep="\t")


df.set_index("Codon", inplace=True)
df_numeric = df.select_dtypes(include=["number"])


plt.figure(figsize=(12, 12))
sns.heatmap(df_numeric, cmap="coolwarm", annot=True, fmt=".2f", cbar_kws={"label": "RSCU Value"})


plt.title("Heatmap of RSCU Values Across Species")
plt.xlabel("Species")
plt.ylabel("Codon")
plt.tight_layout()


plt.savefig("rscu_heatmap.png", dpi=300)
plt.close()
