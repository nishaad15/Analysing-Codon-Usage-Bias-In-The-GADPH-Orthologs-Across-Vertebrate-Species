import pandas as pd

df = pd.read_csv("codon_usage_table.tsv", sep="\t")
df.columns = df.columns.str.strip()
df = df.rename(columns={df.columns[0]: "Codon"})  

long_df = df.melt(id_vars="Codon", var_name="Species", value_name="Count")

codon_to_aa = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser", "CCU": "Pro", "CCC": "Pro",
    "CCA": "Pro", "CCG": "Pro", "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala",
    "GCA": "Ala", "GCG": "Ala", "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop", "CAU": "His",
    "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys", "GAU": "Asp",
    "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys", "UGG": "Trp", "CGU": "Arg", "CGC": "Arg",
    "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg", "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}

long_df["AminoAcid"] = long_df["Codon"].map(codon_to_aa)


def compute_rscu(group):
    total = group["Count"].sum()
    n = len(group)
    group["RSCU"] = (group["Count"] * n) / total if total > 0 else 0
    return group

rscu_long = long_df.groupby(["Species", "AminoAcid"]).apply(compute_rscu).reset_index(drop=True)


rscu_wide = rscu_long.pivot_table(index=["Codon", "AminoAcid"], columns="Species", values="RSCU").reset_index()


rscu_wide = rscu_wide.sort_values("AminoAcid")


rscu_wide.to_csv("rscu_table_wide.tsv", sep="\t", index=False)
print("Saved RSCU table with species as columns to rscu_table_wide.tsv")
