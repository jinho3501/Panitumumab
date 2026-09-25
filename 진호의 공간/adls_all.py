import pandas as pd

adls = pd.read_csv('../csv_file/ADLS_PDS2019.csv')
final_kras = pd.read_csv('../csv_file/final_kras.csv')

target = adls[adls['LSCAT']=='Target lesion'].copy()

reader_rank = {'Radiologist 1': 1,'Radiologist 2':2,'Oncologist':3}
target['_rank'] = target['LSREADER'].map(reader_rank)

target = (target.sort_values(['SUBJID', 'VISITDY', '_rank'])
                .groupby(['SUBJID', 'VISIT', 'VISITDY'], as_index=False)
                .first())

base = target[target['VISIT'] == 'Screening'][['SUBJID', 'LSSLD']]
base = base.rename(columns={'LSSLD': 'BASE'})

sld = target[target['VISIT'] != 'Screening'][['SUBJID', 'VISIT', 'VISITDY', 'LSSLD']]
sld = sld.rename(columns={'LSSLD': 'SLD'})

final_adls_all = sld.merge(base, on='SUBJID', how='left')
final_adls_all['CHG'] = final_adls_all['SLD'] - final_adls_all['BASE']

final_adls_all = final_adls_all.merge(
    final_kras[['SUBJID', 'TRT', 'TRT_bin', 'KRAS_bin']],
    on='SUBJID', how='inner'
)

final_adls_all = final_adls_all[final_adls_all['BASE'].notna()].copy()

print("환자수:", final_adls_all['SUBJID'].nunique())
print("행수:", len(final_adls_all))

final_adls_all.to_csv('../csv_file/adls_all_visits.csv', index=False)