import pandas as pd
import seaborn as sns
import json
import matplotlib.pyplot as plt


def outlier_extraction():

    assemblers = ["Null", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM",
                  "Raven"]
    assemblers_nr = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    assemblers_nr = [str(x) for x in assemblers_nr]

    df = []

    # Table for # of Deletions and Insertions
    for i in range(1, 13):
        counter_deletions = 0
        counter_insertions = 0
        file_path = (f'/home/noah/Research_Project_Neher/data/assembler_comparision/'
                     f'cigar_extraction_files/{i}.json')
        data_insertions_outliers = []
        with open(file_path, 'r') as file:
            json_data = json.load(file)
            value_deletions = json_data.get('Deletions', None)
            value_insertions = json_data.get('Insertions', None)

            for v_i in value_insertions:
                if v_i > 10:
                    data_insertions_outliers.append([assemblers_nr[i], assemblers[i], v_i])

            df_i = pd.DataFrame(data_insertions_outliers, columns=["Assembler Nr.", "Assembler", "Outlier Size"])

            #print("Insertions")
            #print("JSON-File" + str(i))
            #print(df_i)
            df.append(df_i)
    df = pd.concat(df)
    return df



def plot_outliers(df):
    sns.set_style("darkgrid")
    plt.figure(figsize=(10, 6))
    df = df.sort_values("Assembler")
    palette = sns.color_palette("husl", df['Assembler'].nunique())
    ax = sns.stripplot(x='Assembler Nr.', y='Outlier Size', data=df, hue="Assembler", palette=palette, jitter=True)

    ax.set_title('Outliers within insertion sequences', fontsize=16, fontweight='bold')
    ax.set_xlabel('Assembly number', fontsize=14, fontweight='bold')
    ax.set_ylabel('Length of outlier', fontsize=14, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(title='Assemblers', title_fontsize='13', loc='upper left', fontsize='12')
    plt.tight_layout()
    plt.savefig(f"figs/outliers/outlier.png", dpi=300)
    plt.close()


def main():
    df = outlier_extraction()
    plot_outliers(df)


if __name__ == "__main__":
    main()

