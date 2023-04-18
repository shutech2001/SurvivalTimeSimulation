%let dir = C:\Users\Shimizu\国立研究開発法人 国立国際医療研究センター\biostat share - General\Mpox_RCT;

ods html close;

/*
シミュレーションデータ生成用のマクロ

DATA : シミュレーション用に作成するデータ名
SIM : シミュレーション回数
N : 例数（両群合わせたもの．半分ずつ各群に分けられる．）
CONTROL_PROP : 対照群の発症確率
VE : vaccine efficacy
MONTH : 観察期間
*/
%macro DataGenerate(DATA=, SIM=, N=, CONTROL_PROP=, VE=, MONTH=);
data &DATA.;
  call streaminit(42);
  do SIM = 1 to &SIM.;
    do ID = 1 to &N.;
      CONTROL_PROP = &CONTROL_PROP.; VE = &VE.;
	  /* 対照群 */
      if ID <= &N./2 then do;
	    GRP = 0; P = CONTROL_PROP; end;
	  else do;
	    GRP = 1; P = (1 - VE) / 100; end;
      /* 累積確率 = p より 1/lambdaを計算 */
      TIME = rand("exponential", -(356.25/12) * &MONTH./log(1-P));
	  if TIME >= (356.25/12) * &MONTH. then do;
        TIME = (356.25/12) * &MONTH.;
	    CENSOR = 1; EVENT = 0;
	  end;
	  else do;
	    CENSOR = 0; EVENT = 1;
	  end;
    output;end;end;
run;
%mend;

/*
P値を計算するためのマクロ
リスク差・ログランク検定・ポアソン回帰の3つの解析を行う

DATA : 解析対象となるデータセット名
*/
%macro CalcPvalue(DATA=);
  /* リスク差 */
  proc freq data = &DATA.;
    by sim;
    table grp*event / riskdiff(equal column=1 var=sample) chisq nocol nopercent;
    ods output PdiffTest=chisq(where=(Name1="P2_RDIF1"));
  run;
  /* ログランク */
  proc lifetest data = &DATA.;
    by sim;
    time time*censor(1);
    strata / group = grp;
    ods output HomTests=logrank(where=(Test="ログランク"));
  run;
  /* ポアソン回帰 */
  proc genmod data=&DATA.;
    by sim;
    class grp;
    model event(ref="0") = grp / dist = poisson link = log offset = TIME;
    ods output ParameterEstimates=poisson(where=(Parameter="GRP" and Level1 = "0"));
  run;
%mend;

/*
alpha error, beta errorを計算するためのマクロ
出力データはエラーの値が1つだけ格納されたデータとなる

DATA : 解析結果が格納されたデータ
OUT : エラー率が格納される出力データの名前
SIM : シミュレーション回数（エラーを算出するために必要）
PVAL : 解析結果データでP値の列の名前
ALPHA : 有意水準
H0 : H0が正しいか否か（if H0 is TRUE then calculate alpha error, otherwise calculate beta error (power)）
*/
%macro calcErr(DATA=, OUT=, SIM=, PVAL=, ALPHA=, H0=);
  data &DATA.;
    set &DATA.;
	if &PVAL. < &ALPHA. then rej_flg=1; else rej_flg=0; run;
  /* エラーフラグ立てる */
  data flg;
    set &DATA.;
	retain rej 0; rej+rej_flg; *rejectを数え上げる (for alpha error);
    retain acc 0; acc+abs(1-rej_flg); *acceptを数え上げる (for beta error);
	%if &H0.=TRUE %then %do; alpha_error = rej/&sim.; %end;
	%else %do; beta_error = acc/&sim.; power = 1 - beta_error; %end;
  run;
  proc sort data=flg;
    %if &H0. = TRUE %then %do; by descending alpha_error; %end;
    %else %do; by descending beta_error; %end;
  run;
  data &OUT.; set flg(obs=1); 
  %if &H0. = TRUE %then %do; keep alpha_error; %end;
  %else %do; keep power; %end;
  run;
%mend;

/*
DATA : データ
FILEPATH : ファイルのパス
SHEETNAME : シートの名前
*/
%macro out_xlsx(DATA=, FILEPATH=, SHEETNAME=);
  proc export data = &data. outfile = &filepath. dbms = xlsx replace;
    sheet = &sheetname.;
    label;
  run;
%mend;

/*
出力形式を整えたテーブルを出力するマクロ

DATA : シミュレーション用データの名前
SIM : シミュレーション回数
N : 例数（両群合わせたもの）
CONTROL_PROP : 対照群の発症確率
VE : vaccine efficacy
MONTH : 観察期間（月）
H0 : H0が正しいか否か
ALPHA : 有意水準
SET : すでにあるデータセットに縦結合するか否か（2つ目のシナリオ以降はもともとあるテーブルに縦結合していく）
*/
%macro createErrorTab(DATA=, SIM=, N=, CONTROL_PROP=, VE=, MONTH=, H0=, ALPHA=, SET=);
  %DataGenerate(DATA=&DATA., SIM=&SIM., N=&N., CONTROL_PROP=&CONTROL_PROP., VE=&VE., MONTH=&MONTH.);
  %CalcPvalue(DATA=&DATA.);
  %calcErr(DATA=chisq, OUT=out1, SIM=&SIM., PVAL=nValue1, ALPHA=&ALPHA., H0=&H0.);
  %calcErr(DATA=logrank, OUT=out2, SIM=&SIM., PVAL=ProbChiSq, ALPHA=&ALPHA., H0=&H0.);
  %calcErr(DATA=poisson, OUT=out3, SIM=&SIM., PVAL=ProbChiSq, ALPHA=&ALPHA., H0=&H0.);
  data partTab;
    format MONTH CONTROL_PROP VE prop_of_rej_chisq prop_of_rej_logrank prop_of_rej_poisson;
    merge out1(rename=(power=prop_of_rej_chisq)) out2(rename=(power=prop_of_rej_logrank)) out3(rename=(power=prop_of_rej_poisson));
    MONTH = &MONTH.; CONTROL_PROP = &CONTROL_PROP.; VE = &VE.;
  run;
  /* 縦結合 */
  %if &SET. = TRUE %then %do;
    data outTab; set outTab partTab; run;
  %end;
  %if &SET. = FALSE %then %do;
    data outTab; set partTab; run;
  %end;
%mend;
/*
2週間毎にイベントが何件発生しているか分かるようなアウトプット
総数(all_sum)とgroup毎(grp_sum) 

DATA : シミュレーション用のデータセット名
MONTH : 観察期間
CONTROL_PROP : 対照群の発症確率(シート名に記載するため)
VE : vaccine efficacy（シート名に記載するため）
*/
%macro eventEvery2weeks(DATA=, MONTH=, CONTROL_PROP=, VE=);
  %do i=1 %to %sysfunc(int(((356.25/12) * &MONTH.) / 14 + 1));
    data tmpdata;
	  set &DATA.;
      /* 2週間ずつ切り出す */
	  where TIME>=14*(&i.-1) & TIME<=14*&i.;
	proc means sum;
	  by sim;
	  class grp;
	  var event;
	  ods output Summary=sum&i.;
	run;
	/* split data */
	data grp0_sum&i.;
	  set sum&i.; where grp=0; keep sim event_sum;
	data grp1_sum&i.;
	  set sum&i.; where grp=1; keep sim event_sum;
	data sum&i.;
	  merge grp0_sum&i.(rename=(event_sum=event_grp0_sum&i.)) grp1_sum&i.(rename=(event_sum=event_grp1_sum&i.));
	  by sim;
	  event_sum&i. = event_grp0_sum&i. + event_grp1_sum&i.;
	  label event_sum&i. = "全イベント数：期間&i." event_grp0_sum&i. = "対照群イベント数：期間&i." event_grp1_sum&i. = "試験群イベント数：期間&i.";
	run;
	/* merged data */
    /* 2期間目以降は横結合 */
    %if &i.=1 %then %do;
      data all_sum;
        set sum&i.(keep=sim event_sum&i.);
      run;
	  data grp_sum;
	    set sum&i.(keep=sim event_grp0_sum&i. event_grp1_sum&i.);
	  run;
    %end;
    %else %do;
      data all_sum;
	    merge all_sum sum&i.(keep=sim event_sum&i.); by sim;
	  run;
	  data grp_sum;
	    merge grp_sum sum&i.(keep=sim event_grp0_sum&i. event_grp1_sum&i.); by sim;
	  run;
    %end;
  %end;
  %out_xlsx(DATA=all_sum, FILEPATH="&dir.\玉野\eventEvery2weeks.xlsx", SHEETNAME="全例 CONTROL_PROP = &CONTROL_PROP. VE = &VE. MONTH = &MONTH.");
  %out_xlsx(DATA=grp_sum, FILEPATH="&dir.\玉野\eventEvery2weeks.xlsx", SHEETNAME="グループ別 CONTROL_PROP = &CONTROL_PROP. VE = &VE. MONTH = &MONTH.");
%mend;

/*
CONTROL_PROP, VEの値を動的に変えてシミュレーションを実行するためのマクロ

%doマクロの仕様上，ステップ幅はint型でないといけないため使い方要注意

e.g., CONTROL_PROPを0.005 ~ 0.010まで0.001ずつ変えてシミュレーションを行いたいとき
  ==> CTRL_MIN = 5, CTRL_MAX = 10, CTRL_STEP = 1, CTRL_E = 1000
         i.e., 全て整数値になるようにCTRL_E倍する（VEの場合も同様）

<< シミュレーション回数等を変更する場合は，マクロ内のSIM等を変更 >>

CTRL_MIN : CONTROL_PROPの最小値 / CTRL_E
CTRL_MAX : CONTROL_PROPの最大値 / CTRL_E
CTRL_STEP : CONTROL_PROPのステップ幅 / CTRL_E
CTRL_E : 
VE_MIN : VEの最小値 / VE_E
VE_MAX : VEの最大値 / VE_E
VE_STEP : VEのステップ幅 / VE_E
VE_E : 
MONTH : 観察期間（シート名としている）
*/
%macro Validate_CTRL_VE(CTRL_MIN=, CTRL_MAX=, CTRL_STEP=, CTRL_E=, VE_MIN=, VE_MAX=, VE_STEP=, VE_E=, MONTH=);
  %do p_ = &CTRL_MIN. %to &CTRL_MAX. %by &CTRL_STEP.;
    %let p = %sysevalf(&p_./&CTRL_E.);
    %do ve_ = &VE_MIN. %to &VE_MAX. %by &VE_STEP.;
      %let ve = %sysevalf(&ve_./&VE_E.);
      %if &p_. = &CTRL_MIN. and &ve_. = &VE_MIN. %then %createErrorTab(DATA=data, SIM=1000, N=5000, CONTROL_PROP=&p., VE=&ve., MONTH=&MONTH., H0=FALSE, ALPHA=.05, SET=FALSE);
      %else %createErrorTab(DATA=data, SIM=1000, N=5000, CONTROL_PROP=&p., VE=&ve., MONTH=&MONTH., H0=FALSE, ALPHA=.05, SET=TRUE);
      /* 2週毎のイベント数を出すときは下のコメントアウトも実行（実行時間かかる） */
      *%eventEvery2weeks(DATA=data, MONTH=&MONTH., CONTROL_PROP=&p., VE=&ve.);
    %end;
  %end;
  %out_xlsx(data=outTab, filepath="&dir.\玉野\error.xlsx", sheetname="MONTH=&MONTH.");
  /* outTabを削除して意識せずに既にあるデータに縦結合することを防ぐ */
  proc delete data=outTab; run;
%mend;

/* CONTROL_PROP = 0.01, VE = 0.5~0.8の例 */
%Validate_CTRL_VE(CTRL_MIN=1, CTRL_MAX=1, CTRL_STEP=1, CTRL_E=100, VE_MIN=5, VE_MAX=8, VE_STEP=1, VE_E=10, MONTH=2);
