%LET path = /home/student/project_db/csv;

proc import datafile="&path/ADSL_PDS2019.csv"
			out = work.adsl
			dbms= csv replace;
			getnames=yes;
run;

proc import datafile="&path/BIOMARK_PDS2019.csv"
			out= work.biomark
			dbms= csv replace;
			getnames=yes;
run;

proc sort data=work.biomark
	out=work.bio_dedup nodupkey;
	by subjid bmmtr1;
run;

proc freq data=work.bio_dedup;
	tables BMMTR1;
run;

/* CREATE analysis dataset */

proc sort data=work.adsl;
	by subjid;
run;

proc sort data=work.bio_dedup;
	by subjid;
run;

data work.surv;
	merge work.adsl(in=a) work.bio_dedup(keep=subjid bmmtr1);
	by subjid;
	if a;
	* 치료군: 0=BSC, 1=Panitumumab;
	TRTN = (TRT= 'panit. plus best supportive care');

	* KRAS : 0=wild-type, 1=Mutant, Failure은 결측 -> 모델에서 자동 제외 ;
	if BMMTR1 = 'Wild-type' then KRASN = 0;
	else if BMMTR1 = 'Mutant' then KRASN = 1;
	else KRASN = .; *failure 29명;

	KRASEVAL = (KRASN ne .);

	* PFS 시간 0 처리 ;
	PFS0FL = (PFSDYCR <= 0);

	*그대로 두고 제외되게 함;
	PFSTIME = PFSDYCR;

/* 1일로 대체하여 표본에 남김 (주석제거)
if PFSDYCR <= 0 then PFSTIME =1;
else PFSTIME = PFSDYCR ;
*/

	OSTIME = DTHDYX;

label   TRTN='치료군(0=BSC, 1=Panit)'
		KRASN='KRAS(0=WT, 1=Mut)'
		PFSTIME='PFS 기간(일)'
		OSTIME='OS 기간(일)';
run;

proc freq data=work.surv;
	tables BMMTR1 * KRASN /
	missing nopercent nocol;
run;

* 1행=1환자 유지됐는지 확인 (370/370) ;
proc sql;
	select count(*) as n_rows, count(distinct SUBJID) as n_subj from work.surv;
quit;

* KRAS 매핑 제대로 됐는지 (대소문자)
기대: Wild-type 195 / Mutant 146 / 합계 341 ;
proc freq data=work.surv;
	tables BMMTR1 * KRASN / missing nopercent nocol;
run;

* 이벤트 플래그 방향 (결정 7검증)
기대: PFSCR=1이 345명(93%), DTHX=1이 335명(90.5%)
만약 1이 소수라면 1=중도절단.;
proc freq data=work.surv;
	tables PFSCR DTHX PFS0FL KRASEVAL;

proc means data=work.surv n sum median maxdec=1;
	class TRTN KRASN;
	var PFSCR PFSTIME DTHX OSTIME;
	where KRASEVAL;
run;

/* Kaplan-Meier (SAP 6.1) 
Cox보다 먼저 그려서 곡선이 실제로 어떻게 생겼는지 확인 */

* (1) 전체집단 : 치료군 비교;
proc lifetest data=work.surv plots=survival(atrisk=0 to 400 by 50);
	time PFSTIME*PFSCR(0); *괄호 안 0=중도절단 값;
	strata TRTN;
	title 'Figure 1a - PFS 전체 집단';
run;

proc lifetest data=work.surv plots=survival (atrisk=0 to 800 by 100);
	time OSTIME*DTHX(0);
	strata TRTN;
	title 'Figure 2a - OS 전체 집단';
run;

* (2) KRAS subgroup 별 : Wild-type에서는 두 곡선이 벌어지고
Mutant에서는 겹쳐야 상호작용 가설과 일치 ;
proc sort data=work.surv 
	out=work.surv_k;
	by KRASN;
	where KRASEVAL;
run;

proc lifetest data=work.surv_k plots=survival (atrisk=0 to 400 by 50);
	time PFSTIME*PFSCR(0);
	strata TRTN;
	by KRASN;
	title 'Figure 1b = PFS, KRAS subgroup 별';
run;

proc lifetest data=work.surv_k
	plots=survival(atrisk=0 to 800 by 100);
	time OSTIME*DTHX(0);
	strata TRTN;
	by KRASN;
	title 'Figure 2b - os, KRAS subgroup 별';
run;
title;

/* step 6. COX 비례위험모형 (SAP 6.1 두번째)
- ties=efron : 이 실험은 영상 평가가 week 8,12,16,24에만 예정되어 있어
생존시간이 특정 값에 뭉쳐있다. (동점 다수)
동점이 많을 때 Efron 방식이 기본값(Breslow)보다 정확하다.
- hazardratio ... at(KRASN=ALL) : KRAS 수준별 치료 HR을 각각 출력. 6.3 forest plot에 그대로 들어갈 숫자 */

* PFS 1차 평가변수 ;
proc phreg data=work.surv plots(overlay)=survival;
	class TRTN(ref='0') KRASN(ref='0') / param=ref;
	model PFSTIME*PFSCR(0) = TRTN KRASN TRTN*KRASN / rl ties=efron;
	hazardratio 'PFS: KRAS 수준별 치료효과' TRTN / at(KRASN=ALL);
	title 'Table3- Cox PFS (TRT x KRAS)';
run;

* OS 2차 평가변수 ;
proc phreg data=work.surv;
	class TRTN(ref='0') KRASN(ref='0') / param=ref;
	model OSTIME*DTHX(0) = TRTN KRASN TRTN*KRASN / rl ties=efron;
	hazardratio 'OS: KRAS 수준별 치료효과' TRTN / at(KRASN=ALL);
	title 'Table3 - Cox OS (TRT x KRAS) / crossover 영향 해석 필요';
run;

/* step 7. 비례위험 가정 검정
- Cox의 전제는 "두 군의 위험 비가 시간이 지나도 일정하다"
- ASSESS ph / resample : 시뮬레이션으로 p-value를 준다.
p < 0.05 면 가정 위배 신호, 시간구간을 나눠 보거나 RMST 같은 대안 언급 */

proc phreg data=work.surv;
	class TRTN(ref='0')
		  KRASN(ref='0') / param=ref;
	model PFSTIME*PFSCR(0) = TRTN KRASN TRTN*KRASN / ties=efron;
	assess ph / resample seed=20260914;
	title '검증4 - PFS 비례위험 가정';
run;

proc phreg data=work.surv;
	class TRTN(ref='0')
		  KRASN(ref='0') / param=ref;
	model OSTIME*DTHX(0) = TRTN KRASN TRTN*KRASN / ties=efron;
	assess ph / resample seed = 20260914;
	title '검증 4 - OS 비례위험 가정';
run;

/* step 8. 공변량 보정 모형 (선택)
- 원 시험은 ECOG 수행상태와 지역으로 층화 배정했고
원 분석도 이를 보정했다. 현재 SAP에는 공변량이 없으므로
보정모형 함께 제시, 결과가 안정적인지 제시 */

proc phreg data=work.surv;
	class TRTN(ref='0')
		  KRASN(ref='0')
		  B_ECOG SEX / param=ref;
	model PFSTIME*PFSCR(0) = TRTN KRASN TRTN*KRASN AGE B_ECOG SEX / rl ties=efron;
	hazardratio TRTN / at(KRASN=ALL);
	title 'Table 3 보조 - 공변량 보정 PFS';
run;