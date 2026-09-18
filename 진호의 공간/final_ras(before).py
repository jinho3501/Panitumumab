import pandas as pd

bio = pd.read_csv('../csv_file/BIOMARK_FINAL.csv')
final = pd.read_csv('../csv_file/adsl_responder.csv')

final = final.merge(
    bio[['SUBJID','BMMTR1','BMMTR2','BMMTR3','BMMTR4','BMMTR5','BMMTR6']],
    on='SUBJID',how='left'
)
bm_ras = ['BMMTR1','BMMTR2','BMMTR3','BMMTR4','BMMTR5','BMMTR6']

has_failure = (final[bm_ras] == 'Failure').any(axis=1)
has_mutant  = (final[bm_ras] == 'Mutant').any(axis=1)

final['RAS_status'] = 'Wild-type'
final.loc[has_mutant,  'RAS_status'] = 'Mutant'
final.loc[has_failure, 'RAS_status'] = None


final_ras = final[final['RAS_status'].notna()].copy()
final_ras['RAS_bin'] = (final_ras['RAS_status'] == 'Wild-type').astype(int)
final_ras['TRT_bin'] = (final_ras['TRT'].str.contains('panit')).astype(int)

final_ras.to_csv("../csv_file/final_ras(before).csv", index=False)