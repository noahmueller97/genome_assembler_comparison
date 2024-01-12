import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt


input_dir = '/home/noah/Research_Project_Neher/data/assembler_comparision/cigar_files'
assemblers = ["Null", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM",
                  "Raven"]


def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=12,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')
    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)
    for k, cell in mpl_table._cells.items():
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0] % len(row_colors) ])
    return ax
def main():
    data = []

    file_paths = [f'{input_dir}/{i}.paf' for i in range(1, 13)]
    for file_number, file_path in enumerate(file_paths, start=1):
        df = pd.read_csv(file_path, sep='\t', header=None)

        line_q = 0
        line_r = 0

        for idx, row in df.iterrows():
            line_q += (row[3] - row[2])
            line_r += (row[8] - row[7])
        total_q = row[1]
        total_r = row[6]
        unmapped_q = total_q - line_q
        unmapped_r = total_r - line_r

        data.append({"Assembly Number": file_number,
                     "Assembler": assemblers[file_number],
                     "Unmapped regions on the query sequence": unmapped_q,
                     "Unmapped regions on the reference sequence": unmapped_r})

    df_summary = pd.DataFrame(data)

    print(df_summary)
    ax = render_mpl_table(df_summary, header_columns=0, col_width=6.0)
    plt.savefig('/home/noah/Research_Project_Neher/data/assembler_comparision/unmapped_regions/table.png',
                bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()
