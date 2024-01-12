import json
import os
import re
from itertools import groupby
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

input_dir = '/home/noah/Research_Project_Neher/data/assembler_comparision/cigar_files'
output_dir = '/home/noah/Research_Project_Neher/data/assembler_comparision/cigar_extraction_files'
snps = []
data_mutations = []


def parse_cigar(cig):
    cig_iter = groupby(cig, lambda c: c.isdigit())
    for g, n in cig_iter:
        yield int("".join(n)), "".join(next(cig_iter)[1])

def main():
    file_paths = [f'{input_dir}/{i}.paf' for i in range(1, 13)]
    for file_path in file_paths:
        df = pd.read_csv(file_path, sep='\t', header=None)
        matches = []
        deletions = []
        insertions = []
        number_of_matching_bases = []
        gap_compressed = []

        primary = True
        for idx, row in df.iterrows():
            if row[16] != "tp:A:P":
                primary = False
                break

        if primary:
            for cigar_data in df[22]:
                if cigar_data.startswith('cg:Z:'):
                    cigar_string = cigar_data.split('cg:Z:')[1]
                    try:
                        for n, c in parse_cigar(cigar_string):
                            if c == "M":
                                matches.append(n)
                            elif c == "I":
                                insertions.append(n)
                            elif c == "D":
                                deletions.append(n)

                    except ValueError as e:
                        print(f"Error processing CIGAR string: {cigar_string}. Error: {e}")
            for idx, row in df.iterrows():
                pattern = r'\d+'
                pattern2 = r'-?\d+\.\d+|-?\d+'
                tmp_nof = re.findall(pattern, str(row[9]))
                tmp_gap = re.findall(pattern2, row[20])
                number_of_matching_bases.append(tmp_nof)
                gap_compressed.append(tmp_gap)
        else:
            matches = []
            insertions = []
            deletions = []

        res = {'Matches': matches,
               'Insertions': insertions,
               'Deletions': deletions,
               'Number of matching bases': number_of_matching_bases,
               'Number of compressed Gaps': gap_compressed}

        file_name = os.path.basename(file_path).split('.')[0]
        output_file = f'{output_dir}/{file_name}.json'
        with open(output_file, "w") as f:
            json.dump(res, f)

        nomb0 = res["Number of matching bases"]
        nomb1 = int(nomb0[0][0])
        nomb2 = int(nomb0[1][0])
        nomb = nomb1 + nomb2


        matches_sum = sum(res["Matches"])
        number_of_matches = len(matches)

        snps.append((matches_sum - nomb)/matches_sum)

    print(snps)
    assemblers = ["Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM", "Raven"]
    for i in range(12):
        data_mutations.append([assemblers[i], snps[i]])
    df_mutations = pd.DataFrame(data_mutations, columns=["Assembler", "SNPS"])

    sns.set_style("darkgrid")
    plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
    palette = sns.color_palette("flare", len(df_mutations['Assembler'].unique()))  # Color palette
    ax = sns.stripplot(x='Assembler', y='SNPS', data=df_mutations, hue="Assembler", palette=palette, jitter=True)
    ax.set_title('Mutation Rate', fontsize=16, fontweight='bold')
    ax.set_xlabel('Assembler', fontsize=14, fontweight='bold')  # Optionally, set the label for the x-axis
    ax.set_ylabel("Mismatches density", fontsize=14, fontweight='bold')  # Optionally, set the label for the y-axis
    ax.tick_params(axis='both', which='major', labelsize=12)  # Adjust the size of the axis ticks
    plt.tight_layout()
    plt.savefig(f"figs/mutations/mutation.png", dpi=300)
    plt.close()



if __name__ == "__main__":
    main()
