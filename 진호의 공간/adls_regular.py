import pandas as pd

adls = pd.read_csv('../csv_file/ADLS_PDS2019.csv')
final_kras = pd.read_csv('../csv_file/final_kras.csv')


target = adls[adls['LSCAT']=='Target lesion'].copy()

# ADRSP처럼 우선순위를 Radiologist 1 > Radiologist 2 > Oncologist
reader_rank = {'Radiologist 1': 1,'Radiologist 2':2,'Oncologist':3}
target['_rank'] = target['LSREADER'].map(reader_rank)

target = (target.sort_values(['SUBJID', 'VISITDY', '_rank'])
                .groupby(['SUBJID', 'VISIT', 'VISITDY'], as_index=False)
                .first())

# 각 정규 방문
regular_visits = ['Screening', 'Week 8', 'Week 12', 'Week 16', 'Week 24',
                   'Week 32', 'Week 40', 'Week 48', 'Week 60',
                   'Week 72', 'Week 84', 'Week 96', 'Week 108', 'Week 120', 'Week 132']
target = target[target['VISIT'].isin(regular_visits)].copy()

#기저치(BASE) : Screening 시점 표적 병변 직경 총합
base = target[target['VISIT'] == 'Screening'][['SUBJID', 'LSSLD']]
base = base.rename(columns={'LSSLD': 'BASE'})

#SLD : 각 정규 방문(Week 8, 12, 16, 24 등) 시점의 표적 병변 직경 총합
sld = target[target['VISIT'] != 'Screening'][['SUBJID', 'VISIT', 'VISITDY', 'LSSLD']]
sld = sld.rename(columns={'LSSLD': 'SLD'})

final_adls = sld.merge(base, on='SUBJID', how='left')

#각 방문 시점의 기저치 대비 변화량 (CHG = SLD - BASE)
final_adls['CHG'] = final_adls['SLD'] - final_adls['BASE']

#ADLS에 기저 후(Post-baseline) 측정이 1회 이상 존재하는 KRAS 평가군
final_adls = final_adls.merge(
    final_kras[['SUBJID', 'TRT', 'TRT_bin', 'KRAS_bin']],
    on='SUBJID', how='inner'
)

final_adls= final_adls[final_adls['BASE'].notna()].copy()

final_adls.to_csv('../csv_file/adls_regular_visits.csv', index=False)