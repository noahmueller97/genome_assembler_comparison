import pandas as pd
import seaborn as sns
import json
import matplotlib.pyplot as plt

data_amount_of_indels = []
data_sum_of_indels = []

# Creates a table for # of Deletions and Insertions
for i in range(1, 13):
    counter_deletions = 0
    counter_insertions = 0
    file_path = (f'/home/noah/Research_Project_Neher/data/assembler_comparision/'
                 f'cigar_extraction_files/{i}.json')
    with open(file_path, 'r') as file:
        json_data = json.load(file)
        value_deletions = json_data.get('Deletions', None)
        value_insertions = json_data.get('Insertions', None)
        for v_d in value_deletions:
            counter_deletions += 1
        for v_i in value_insertions:
            counter_insertions += 1
        assemblers = ["Null", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM", "Raven"]
        data_amount_of_indels.append([i, assemblers[i], counter_deletions, counter_insertions])
    df_amount_of_indels = pd.DataFrame(data_amount_of_indels, columns=["Assembly Nr.", "Assembler", "# of Deletions", "# of Insertions"])

# Table for the sum of Deletions and Insertions
for i in range(1, 13):
    counter_deletions = 0
    counter_insertions = 0
    file_path = (f'/home/noah/Research_Project_Neher/data/assembler_comparision/'
                 f'cigar_extraction_files/{i}.json')
    with open(file_path, 'r') as file:
        json_data = json.load(file)
        value_deletions = json_data.get('Deletions', None)
        value_insertions = json_data.get('Insertions', None)
        for v_d in value_deletions:
            counter_deletions += v_d
        for v_i in value_insertions:
            counter_insertions += v_i
        assemblers = ["Null", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM", "Raven", "Flye", "MM", "Raven"]
        data_sum_of_indels.append([i, assemblers[i], counter_deletions, counter_insertions])
    df_sum_of_indels = pd.DataFrame(data_sum_of_indels, columns=["Assembly Nr.", "Assembler", "Sum of Deletions", "Sum of Insertions"])
print(df_sum_of_indels)
print(df_amount_of_indels)



# Plot of assembler comparison in respect to the number of deletions
sns.set_style("darkgrid")
plt.figure(figsize=(8, 5))
ax = sns.barplot(data=df_amount_of_indels, x='Assembler', y='# of Deletions', palette="pastel")
ax.set_title('Comparison of Assemblers', fontsize=16, fontweight='bold')
ax.set_xlabel('Assembler', fontsize=14, fontweight='bold')
ax.set_ylabel('# of Deletions', fontsize=14, fontweight='bold')

for p in ax.patches:
    height = p.get_height()
    ax.annotate(f'{height}', 
                (p.get_x() + p.get_width() / 2., height), 
                ha='center', va='center', 
                xytext=(-30, 10),  # Adjust this value to move the text higher
                textcoords='offset points',
                fontsize=12, 
                color='black',
                zorder=3)  # Set z-order to draw text on top


plt.savefig(f"figs/assembler_comparison_plots/amount_of_deletions.png")
plt.close()



# Plot of assembler comparison in respect to the number of insertions
sns.set_style("darkgrid")
plt.figure(figsize=(8, 5))
ax = sns.barplot(data=df_amount_of_indels, x='Assembler', y='# of Insertions', palette="pastel")
ax.set_title('Comparison of Assemblers', fontsize=16, fontweight='bold')
ax.set_xlabel('Assembler', fontsize=14, fontweight='bold')
ax.set_ylabel('# of Insertions', fontsize=14, fontweight='bold')

for p in ax.patches:
    height = p.get_height()
    ax.annotate(f'{height}', 
                (p.get_x() + p.get_width() / 2., height), 
                ha='center', va='center', 
                xytext=(-30, 10),  # Adjust this value to move the text higher
                textcoords='offset points',
                fontsize=12, 
                color='black',
                zorder=3)  # Set z-order to draw text on top


plt.savefig(f"figs/assembler_comparison_plots/amount_of_insertions.png")
plt.close()


# Plot of assembler comparison in respect to the sum of deletions
sns.set_style("darkgrid")
plt.figure(figsize=(8, 5))
ax = sns.barplot(data=df_sum_of_indels, x='Assembler', y='Sum of Deletions', palette="pastel")
ax.set_title('Comparison of Assemblers', fontsize=16, fontweight='bold')
ax.set_xlabel('Assembler', fontsize=14, fontweight='bold')
ax.set_ylabel('Total Length Of Deleted Sequence', fontsize=14, fontweight='bold')

for p in ax.patches:
    height = p.get_height()
    ax.annotate(f'{height}', 
                (p.get_x() + p.get_width() / 2., height), 
                ha='center', va='center', 
                xytext=(-30, 10),  # Adjust this value to move the text higher
                textcoords='offset points',
                fontsize=12, 
                color='black',
                zorder=3)  # Set z-order to draw text on top


plt.savefig(f"figs/assembler_comparison_plots/sum_of_deletions.png")
plt.close()

# Plot of assembler comparison in respect to the sum of insertions
sns.set_style("darkgrid")
plt.figure(figsize=(8, 5))
ax = sns.barplot(data=df_sum_of_indels, x='Assembler', y='Sum of Insertions', palette="pastel")
ax.set_title('Comparison of Assemblers', fontsize=16, fontweight='bold')
ax.set_xlabel('Assembler', fontsize=14, fontweight='bold')
ax.set_ylabel('Total Length Of Inserted Sequence', fontsize=14, fontweight='bold')

for p in ax.patches:
    height = p.get_height()
    ax.annotate(f'{height}', 
                (p.get_x() + p.get_width() / 2., height), 
                ha='center', va='center', 
                xytext=(-30, 10),  # Adjust this value to move the text higher
                textcoords='offset points',
                fontsize=12, 
                color='black',
                zorder=3)  # Set z-order to draw text on top


plt.savefig(f"figs/assembler_comparison_plots/sum_of_insertions.png")
plt.close()