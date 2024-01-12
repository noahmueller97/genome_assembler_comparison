import pandas as pd
import re
numbers_of_snps = []
number_of_indels = []

assembler1 = "Flye"
assembler2 = "Miniasm and Minipolish"
assembler3 = "Raven"

assemblers = [0]*12

for x in range(len(assemblers)):
    if x % 3 == 0:
        assemblers[x] = assembler1
    elif x % 3 == 1:
        assemblers[x] = assembler2
    else:
        assemblers[x] = assembler3
assemblers.append("Pre-Medaka")
print(assemblers)


pattern_snps = r"TotalSNPs\s+(\d+)"
pattern_indels = r"TotalIndels\s+(\d+)"
folder_path = '/home/noah/Research_Project_Neher/data/assembler_comparision/reports/'
for i in range(1, 14):
    file_name = f"{folder_path}{i}.report"
    with open(file_name, 'r') as file:
        file_contents = file.read()
        match_snps = re.search(pattern_snps, file_contents)
        match_indels = re.search(pattern_indels, file_contents)
        if match_snps:
            number_snps = match_snps.group(1)
            numbers_of_snps.append(number_snps)
        else:
            print(f"No match found in file: {file_name}")
        if match_indels:
            number_indels = match_indels.group(1)
            number_of_indels.append(number_indels)
        else:
            print(f"No match found in file: {file_name}")

data = {"Assembly": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "SNPs": numbers_of_snps, "Indels": number_of_indels,
        "Assembler": assemblers}
df = pd.DataFrame(data)
df.to_csv("report.csv")
print(df)

