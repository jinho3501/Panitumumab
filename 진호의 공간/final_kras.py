import pandas as pd

bio = pd.read_csv('../csv_file/BIOMARK_PDS2019.csv')
final = pd.read_csv('../csv_file/adsl_responder.csv')

bio_clean = bio.groupby('SUBJID', as_index=False).last()

final = final.merge(
    bio_clean[['SUBJID','BMMTR1']],
    on='SUBJID',how='left'
)

final_kras = final[final['BMMTR1'] != 'Failure'].copy()
final_kras['KRAS_bin'] = (final_kras['BMMTR1'] == 'Wild-type').astype(int)
final_kras['TRT_bin'] = (final_kras['TRT'].str.contains('panit')).astype(int)

print(final_kras)

final_kras.to_csv("../csv_file/final_kras.csv", index=False)
