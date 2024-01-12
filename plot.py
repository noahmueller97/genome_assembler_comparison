import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def reportPlot():
    df = pd.read_csv("report.csv")

    sns.stripplot(data=df, x="Assembler", y="SNPs")
    plt.savefig("figs/medaka/snps.png")
    plt.close()

    sns.stripplot(data=df, x="Assembler", y="Indels")
    plt.savefig("figs/medaka/indels.png")
    plt.close()


def cigarPlot():
    for file_number in range(1,13):
        file_path = (f'/home/noah/Research_Project_Neher/data/assembler_comparision/'
                     f'cigar_extraction_files/{file_number}.json')
        with open(file_path, 'r') as file:
            data = json.load(file)
        target_key_insertions = "Insertions"
        target_key_deletions = "Deletions"
        target_key_matches = "Matches"

        if target_key_insertions in data:
            values_insertions = data[target_key_insertions]
        else:
            print(f"Key '{target_key_insertions}' not found in the JSON file.")
            values_insertions = []

        if target_key_deletions in data:
            values_deletions = data[target_key_deletions]
        else:
            print(f"Key '{target_key_deletions}' not found in the JSON file.")
            values_deletions = []

        if target_key_matches in data:
            values_matches = data[target_key_matches]
        else:
            print(f"Key '{target_key_matches}' not found in the JSON file.")
            values_matches = []

        if values_insertions:
            sns.set_style("darkgrid")
            df_i = pd.DataFrame({'Values': values_insertions})
            bins = np.logspace(0, 3, 20)
            plt.figure(figsize=(10, 6))
            mask = df_i > 0
            ax = sns.histplot(data=df_i[mask], x='Values',
                              log_scale=True,
                              bins=np.log10(bins),
                              stat='count',
                              color='lightgreen')

            # Customizing the plot
            ax.set_title('Insertions', fontsize=16, fontweight='bold')
            ax.set_xlabel('Length of segments (log scale)', fontsize=14, fontweight='bold')
            ax.set_ylabel('Amount of segments', fontsize=14, fontweight='bold')
            ax.set_xticks([1, 2, 3, 4, 5, 10, 20, 40, 100, 250, 600])
            ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())

            plt.tight_layout()
            plt.savefig(f"figs/cigar/Insertion_{file_number}.png", dpi=300, bbox_inches='tight')
            plt.close()

        else:
            print("No data to plot.")


        if values_deletions:
            sns.set_style("darkgrid")
            df_v = pd.DataFrame({'Values': values_deletions})
            plt.figure(figsize=(10, 6))
            ax = sns.histplot(data=df_v, x='Values', bins=range(1, 7), discrete=True,
                              color='skyblue')

            # Customizing the plot
            ax.set_title('Deletions', fontsize=16, fontweight='bold')
            ax.set_xlabel('Length of segments', fontsize=14, fontweight='bold')
            ax.set_ylabel('Amount of segments', fontsize=14, fontweight='bold')
            plt.xticks(range(1, 7), fontsize=12)
            plt.yticks(fontsize=12)

            for p in ax.patches:
                ax.annotate(f'{int(p.get_height())}',
                            (p.get_x() + p.get_width() / 2., p.get_height()),
                            ha='center', va='center',
                            xytext=(0, 5),
                            textcoords='offset points',
                            fontsize=12,
                            color='black')

            plt.savefig(f"figs/cigar/Deletions_{file_number}.png", dpi=300, bbox_inches='tight')
            plt.close()

        else:
            print("No data to plot.")

        if values_matches:
            sns.set_style("darkgrid")
            df_m = pd.DataFrame({'Values': values_matches})
            plt.figure(figsize=(10, 6))

            start = np.log10(min(values_matches))
            stop = np.log10(max(values_matches))
            bins = np.logspace(start, stop, 50)

            ax = sns.histplot(data=df_m, x='Values', bins=bins, color='PaleVioletRed')

            plt.xscale('log')
            plt.xticks([100, 1000, 10000, 100000], ['100', '1k', '10k', '100k'])

            ax.set_title('Matches (Log Scale)', fontsize=16, fontweight='bold')
            ax.set_xlabel('Length of segments', fontsize=14, fontweight='bold')
            ax.set_ylabel('Amount of segments', fontsize=14, fontweight='bold')

            plt.tight_layout()

            plt.savefig(f"figs/cigar/Matches_{file_number}.png", dpi=300, bbox_inches='tight')
            plt.close()

        else:
            print("No data to plot.")


def main():
    cigarPlot()
    reportPlot()

if __name__ == "__main__":
    main()
