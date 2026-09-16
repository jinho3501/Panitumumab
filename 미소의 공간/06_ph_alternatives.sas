
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

/* 공변량을 보고 ph alternatives (대응분석) */
* A. 구간별 위험비 (분할시점-> 112일 = week16);

%let cut = 112;   /* week 16 */

/*--------------------------------------------------------------
  A-1. 한 환자를 두 줄로 쪼개기 (counting process 형식)

  원래 데이터는 환자 1명 = 1행 이다. 이걸 구간마다 한 행씩으로
  쪼개서, 같은 환자가 구간 1과 구간 2에 각각 등장하게 만든다.
"100일까지의 위험"과 "100일 이후의 위험"을
  따로 추정할 수 있다.
--------------------------------------------------------------*/

data work.split;
  set work.surv;
  where KRASEVAL = 1 and PFSTIME > 0;   /* PFS시간 0인 15명은 구간을 만들 수 없어 제외 */

  /*--- 구간 1 : 0 ~ &cut 일 ---*/
  tstart = 0;
  tstop  = min(PFSTIME, &cut);              /* cut 을 넘겼으면 cut 에서 끊는다 */
  event  = (PFSCR = 1 and PFSTIME <= &cut); /* cut 안에서 일어난 사건만 1 */
  period = 1;
  output;

  /*--- 구간 2 : &cut 일 이후 (cut 을 넘긴 환자만) ---*/
  if PFSTIME > &cut then do;
    tstart = &cut;                          /* 이 환자는 100일부터 관찰 시작 */
    tstop  = PFSTIME;
    event  = (PFSCR = 1);
    period = 2;
    output;
  end;

  label period='관찰 구간' tstart='구간 시작일' tstop='구간 종료일';
run;

proc freq data=work.split;
  tables period * event / nopercent nocol;
  title "A-1 검증 — 구간별 행 수와 사건 수 (분할 시점 &cut 일)";
run;


/*--------------------------------------------------------------
  A-2. 구간별로 같은 모형을 따로 적합

  by period;  →  구간 1 결과, 구간 2 결과가 각각 출력된다.

--------------------------------------------------------------*/

proc sort data=work.split; by period; run;

proc phreg data=work.split;
  by period;
  class TRTN(ref='0') KRASN(ref='0') / param=ref;
  model (tstart, tstop)*event(0) = TRTN KRASN TRTN*KRASN / rl ties=efron;
  hazardratio '구간별 KRAS 수준별 치료효과' TRTN / at(KRASN=ALL);
  title "A-2 구간별 Cox — period 1: 0~&cut 일 / period 2: &cut 일 이후";
run;
title;


/*==============================================================
  [B] RMST (Restricted Mean Survival Time, 제한평균생존시간)
==============================================================*/

/*--------------------------------------------------------------
  B-0. RMST 가 무엇인가

  "tau 일까지 평균 며칠을 사건 없이 보냈는가"


  [근거] tau = 168일 = week 24

   (1) 원 시험 프로토콜의 영상 평가 예정 시점이다 (데이터 외부 근거)
   (2) 약 6개월. 6개월 PFS 율은 종양학에서 표준적으로 보고되는 지표다
   (3) 기술적 제약을 만족한다 — tau 는 각 군의 최대 관측시간 중
       작은 값을 넘으면 안 된다. 본 데이터에서 가장 짧은 군도
       400일 이상 관찰되었으므로 168일은 충분히 안쪽이다

  [FDA 지침 참고] 종양학 OS 지침은 비례위험 위반이 예상될 때
   landmark 생존율과 RMST 를 사전에 지정하도록 권고한다. 즉 RMST 는
   임의로 고른 방법이 아니라 지침에 명시된 표준 대응책이다.
--------------------------------------------------------------*/
%let tau = 168;   /* week 24, 약 6개월 */


/*--------------------------------------------------------------
  B-1. KRAS 하위군별 RMST
--------------------------------------------------------------*/

proc lifetest data=work.surv rmst(tau=&tau) plots=none;
  where KRASN = 0;
  time PFSTIME*PFSCR(0);
  strata TRTN;
  title "B-1 RMST — KRAS Wild-type, tau = &tau 일";
run;

proc lifetest data=work.surv rmst(tau=&tau) plots=none;
  where KRASN = 1;
  time PFSTIME*PFSCR(0);
  strata TRTN;
  title "B-1 RMST — KRAS Mutant, tau = &tau 일";
run;


/* 참고 : 전체 KRAS 평가군 */
proc lifetest data=work.surv rmst(tau=&tau) plots=none;
  where KRASEVAL = 1;
  time PFSTIME*PFSCR(0);
  strata TRTN;
  title "B-1 RMST — KRAS 평가군 전체, tau = &tau 일";
run;
title;


/*--------------------------------------------------------------
  B-2. (선택) RMST 회귀모형 — 상호작용을 직접 검정

  PROC RMSTREG 는 RMST 를 결과변수로 하는 회귀모형이다.
  Cox 와 같은 구조로 상호작용을 넣을 수 있고, 계수는 "일(day)"
  단위로 해석된다.
--------------------------------------------------------------*/
proc rmstreg data=work.surv tau=&tau;
  where KRASEVAL = 1;
  class TRTN(ref='0') KRASN(ref='0');
  model PFSTIME*PFSCR(0) = TRTN KRASN TRTN*KRASN / link=linear;
  title "B-2 RMST 회귀 — 상호작용 포함, tau = &tau 일";
run;
title;
