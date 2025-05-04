from Bio import SeqIO
from collections import Counter
import os

def count_codons(seq):
   
    codons = [str(seq[i:i+3]).upper() for i in range(0, len(seq)-2, 3) if len(seq[i:i+3]) == 3]
    return Counter(codons)

def analyze_folder(folder_path, output_file):
    all_results = {}
    
    for filename in os.listdir(folder_path):
        if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".txt"):
            species = os.path.splitext(filename)[0] 
            codon_counter = Counter()
            file_path = os.path.join(folder_path, filename)
            
            for record in SeqIO.parse(file_path, "fasta"):
                codon_counter.update(count_codons(record.seq))
            
            all_results[species] = codon_counter

    
    all_codons = sorted({codon for counters in all_results.values() for codon in counters})
    
    with open(output_file, "w") as out:
        
        out.write("Codon\t" + "\t".join(all_results.keys()) + "\n")
        
        for codon in all_codons:
            row = [codon] + [str(all_results[sp].get(codon, 0)) for sp in all_results]
            out.write("\t".join(row) + "\n")

    print(f"Codon usage written to '{output_file}'.")


folder = "sequences"  
output = "codon_usage_table.tsv"
analyze_folder(folder, output)
