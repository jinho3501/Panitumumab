import pandas as pd

adsl = pd.read_csv("csv_file/ADSL_PDS2019.csv")
adrsp = pd.read_csv("csv_file/ADRSP_PDS2019.csv")

adrsp = adrsp[(adrsp['RSRESP'] != 'Unknown') & (adrsp['VISIT'] != 'Screening')].copy()

# Radiologist 1 으로 기준
reader_rank = {'Radiologist 1': 1,'Radiologist 2':2,'Oncologist':3}
adrsp['_rank'] = adrsp['RSREADER'].map(reader_rank)

adrsp = (adrsp.sort_values(by=['SUBJID','VISITDY','_rank'])
         .groupby(['SUBJID', 'VISITDY'], as_index=False)
          .first())
resp_rank = {
    "Complete response":   1, # 없지만 그냥 넣어놈(국제 표준이라고해서)
    "Partial response":    2,
    "Stable disease":      3,
    "Progressive disease": 4,
    "Unable to evaluate":  5,
}
adrsp['_resp'] = adrsp['RSRESP'].map(resp_rank)

# BOR로 가장 좋았던 판정으로 판정
adrsp = (adrsp.sort_values(['SUBJID', '_resp'])
            .groupby('SUBJID', as_index=False)
            .first()[['SUBJID', 'RSRESP']]
            .rename(columns={'RSRESP': 'BOR'}))

adrsp['RESPONDER'] = adrsp['BOR'].isin(['Complete response', 'Partial response']).astype(int)

final = adsl.merge(adrsp, on='SUBJID', how='left')

final['BOR'] = final['BOR'].fillna('unknown')

# 43명이 평가를 못 받아서 0으로
final['RESPONDER'] = final['RESPONDER'].fillna(0).astype(int)

final.to_csv("../csv_file/adsl_responder.csv", index=False)