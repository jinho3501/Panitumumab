import pandas as pd

final_bio = pd.read_csv("../csv_file/BIOMARK_FINAL.csv")
final_ras = pd.read_csv("../csv_file/final_ras.csv")
before_ras = pd.read_csv("../csv_file/final_ras(before).csv")
final_kras = pd.read_csv("../csv_file/final_kras.csv")
adae = pd.read_csv("../csv_file/ADAE_PDS2019.csv")

print(adae['AESOC'].value_counts())


#
# print(len(final_ras))
# print(len(before_ras))
# print(final_ras['RAS_status'].value_counts())
# print(before_ras['RAS_status'].value_counts())
# print(len(final_kras))
# print(final_kras['BMMTR1'].value_counts())
# print(pd.crosstab(final_ras['TRT_bin'], final_ras['RAS_bin']))

# print(len(final_bio))
# print(final_bio['RAS_status'].value_counts())

