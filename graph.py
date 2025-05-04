import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


df = pd.read_csv("codon_usage_table.tsv", sep="\t")


codon_to_aa = {
    # Phenylalanine
    "UUU": "Phe", "UUC": "Phe",
    # Leucine
    "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    # Isoleucine
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",
    # Methionine (Start)
    "AUG": "Met",
    # Valine
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    # Serine
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser",
    # Proline
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    # Threonine
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    # Alanine
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    # Tyrosine
    "UAU": "Tyr", "UAC": "Tyr",
    # Stop codons
    "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
    # Histidine
    "CAU": "His", "CAC": "His",
    # Glutamine
    "CAA": "Gln", "CAG": "Gln",
    # Asparagine
    "AAU": "Asn", "AAC": "Asn",
    # Lysine
    "AAA": "Lys", "AAG": "Lys",
    # Aspartic Acid
    "GAU": "Asp", "GAC": "Asp",
    # Glutamic Acid
    "GAA": "Glu", "GAG": "Glu",
    # Cysteine
    "UGU": "Cys", "UGC": "Cys",
    # Tryptophan
    "UGG": "Trp",
    # Arginine
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    # Glycine
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


df["AminoAcid"] = df["Codon"].map(codon_to_aa)
df["Codon"] = df["Codon"].str.upper().str.strip()
df["AminoAcid"] = df["Codon"].map(codon_to_aa)
df = df.dropna(subset=["AminoAcid"])


output_dir = "codon_usage_by_aa"
os.makedirs(output_dir, exist_ok=True)

unique_aas = sorted(df["AminoAcid"].unique())

colors = sns.color_palette("Set2", n_colors=len(df.columns[1:-1]))

fig, axes = plt.subplots(5, 4, figsize=(20, 18))
axes = axes.flatten()

for i, aa in enumerate(unique_aas):
    aa_df = df[df["AminoAcid"] == aa]
    ax = axes[i]

   
    aa_df_melt = aa_df.melt(id_vars=["Codon"], value_vars=df.columns[1:-1], 
                            var_name="Species", value_name="Frequency")
    
    sns.barplot(data=aa_df_melt, x="Codon", y="Frequency", hue="Species", ax=ax, palette=colors)
    ax.set_title(f"{aa}", fontsize=14)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.legend().set_visible(False)  

for j in range(len(unique_aas), len(axes)):
    fig.delaxes(axes[j])

#legend
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(labels), fontsize=12)

plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.suptitle("Codon Usage Grouped by Amino Acid Across Species", fontsize=18)
plt.savefig("codon_usage_by_amino_acid.png", dpi=300)
plt.show()
