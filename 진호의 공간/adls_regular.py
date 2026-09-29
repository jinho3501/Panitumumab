import pandas as pd

adls = pd.read_csv('../csv_file/ADLS_PDS2019.csv')
final_kras = pd.read_csv('../csv_file/final_kras.csv')

target = adls[adls['LSCAT']=='Target lesion'].copy()

reader_rank = {'Radiologist 1': 1,'Radiologist 2':2,'Oncologist':3}
target['_rank'] = target['LSREADER'].map(reader_rank)

target = (target.sort_values(['SUBJID', 'VISIT', 'VISITDY', '_rank'])
                .groupby(['SUBJID', 'VISIT'], as_index=False)
                .first())

regular_visits = ['Screening', 'Week 8', 'Week 12', 'Week 16', 'Week 24',
                   'Week 32', 'Week 40', 'Week 48', 'Week 60',
                   'Week 72', 'Week 84', 'Week 96', 'Week 108', 'Week 120', 'Week 132']
target = target[target['VISIT'].isin(regular_visits)].copy()

base = target[target['VISIT'] == 'Screening'][['SUBJID', 'LSSLD']]
base = base.rename(columns={'LSSLD': 'BASE'})

sld = target[target['VISIT'] != 'Screening'][['SUBJID', 'VISIT', 'VISITDY', 'LSSLD']]
sld = sld.rename(columns={'LSSLD': 'SLD'})

final_adls = sld.merge(base, on='SUBJID', how='left')
final_adls['CHG'] = final_adls['SLD'] - final_adls['BASE']

final_adls = final_adls[final_adls['SLD'].notna()].copy()

final_adls = final_adls.merge(
    final_kras[['SUBJID', 'TRT', 'TRT_bin', 'KRAS_bin']],
    on='SUBJID', how='inner'
)

final_adls = final_adls[final_adls['BASE'].notna()].copy()

print("환자수:", final_adls['SUBJID'].nunique())
print("행수:", len(final_adls))


dup = final_adls.groupby(['SUBJID','VISIT']).size()
print("이상 중복 건수:", (dup>1).sum())

final_adls.to_csv('../csv_file/adls_regular.csv', index=False)