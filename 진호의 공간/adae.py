import pandas as pd

adae = pd.read_csv('../csv_file/ADAE_PDS2019.csv')
adsl = pd.read_csv('../csv_file/ADSL_PDS2019.csv')

first = adae.groupby('SUBJID')['AESTDYI'].min().median() # 11.5일 나옴

skin_all = adae[adae['AESOC'] == 'SKIN AND SUBCUTANEOUS TISSUE DISORDERS']

landmark = 15
final_skin_lm=  adsl[adsl['DTHDYX']>=landmark].copy() #6명이 15일 전에 사망

skin_before = skin_all[skin_all['AESTDYI']<=landmark] #9명이 15일 이후에 피부발진 시작
grade_lm = skin_before.groupby('SUBJID')['AESEVCD'].max().reset_index()
grade_lm.columns = ['SUBJID', 'MAX_GRADE_LM']

final_skin_lm = final_skin_lm.merge(grade_lm, on='SUBJID', how='left')
final_skin_lm['MAX_GRADE_LM'] = final_skin_lm['MAX_GRADE_LM'].fillna(0)
final_skin_lm['MAX_GRADE_LM'] = final_skin_lm['MAX_GRADE_LM'].astype(int)
final_skin_lm['SEVERE_RASH_LM'] = (final_skin_lm['MAX_GRADE_LM'] >= 2).astype(int)

# ⑥ 생존시간 재계산 (랜드마크 기준으로)
final_skin_lm['DTHDYX_LM'] = final_skin_lm['DTHDYX'] - landmark

# ⑦ 파니투무맙군만
final_skin_lm = final_skin_lm[final_skin_lm['TRT'].str.contains('panit')]

print(final_skin_lm)
print(len(final_skin_lm))
print(final_skin_lm['SEVERE_RASH_LM'].value_counts())