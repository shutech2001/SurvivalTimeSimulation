%let dir = C:\Users\Shimizu\���������J���@�l �������ۈ�Ì����Z���^�[\biostat share - General\Mpox_RCT;

ods html close;



control_prop VE prop_of_rej_chesq prop_of_rej_logrank prop_of_rej_poisson
sheetname = "control_prop=, VE=, Month="

%macro DataGenerate(DATA=, SIM=, N=, ctrl_rate=, vac_eff=, month=);
data &DATA.;
  call streaminit(42);
  do sim = 1 to &SIM.;
    do id = 1 to &N.;
	  /* �ΏƌQ */
      if id <= &N./2 then do;
	    grp = 0;
        rate = &ctrl_rate.;
	  end;
	  else do;
	    grp = 1;
	    rate = (1 - &vac_eff.) / 100;
	  end;
      /* �ݐϊm�� = rate ��� 1/lambda���v�Z */
      TIME = rand("exponential", -(356.25/12) * &month./log(1-rate));
	  if TIME >= (356.25/12) * &month. then do;
        TIME = (356.25/12) * &month.;
	    censor = 1;
	    event = 0;
	  end;
	  else do;
	    censor = 0;
        event = 1;
	  end;
    output;end;
  end;
run;
%mend;

%macro CalcPvalue(DATA=);
  proc freq data = &data.;
    by sim;
    table grp*event / riskdiff(equal column=1 var=sample) chisq nocol nopercent;
    ods output PdiffTest=chisq(where=(Name1="P2_RDIF1"));
  run;

  proc lifetest data = &data.;
    by sim;
    time time*censor(1);
    strata / group = grp;
    ods output HomTests=logrank(where=(Test="���O�����N"));
  run;

  proc genmod data=&data.;
    by sim;
    class grp;
    model event = grp / dist = poisson link = log offset = TIME;
    ods output ParameterEstimates=poisson(where=(ProbChiSq^=. and Parameter="grp"));
  run;
%mend;

/*
data : �f�[�^
pval : p�l�̗�
alpha : �L�Ӑ���
H0 : �A�����������������ۂ�
*/
%macro calcErr(data=, sim=, pval=, alpha=, H0=);
  data &data.;
    set &data.;
	if &pval. < &alpha. then rej_flg=1; else rej_flg=0; run;
  /* �G���[�t���O���Ă� */
  data flg;
    set &data.;
	retain rej 0;
	rej+rej_flg;
	retain acc 0;
	acc+abs(1-rej_flg);
	%if &H0.=TRUE %then %do; alpha_error = rej/&sim.; keep alpha_error; %end;
	%else %do; beta_error = acc/&sim.; keep beta_error; %end;
  run;
  proc sort data=flg;
    %if &H0. = TRUE %then %do; by descending alpha_error; %end;
    %else %do; by descending beta_error; %end;
  run;
  /* �G���[�o�� */
  data out; set flg(obs=1); run;
  /* proc print data=out; run; */
%mend;
/*
data : �f�[�^
filepath : �t�@�C���̃p�X
sheetname : �V�[�g�̖��O
*/
%macro out_xlsx(data=, filepath=, sheetname=);
  proc export data = &data. outfile = &filepath. dbms = xlsx replace;
    sheet = &sheetname.;
    label;
  run;
%mend;

/* 2�T�Ԃ����ɃC�x���g�������������Ă��邩������悤�ȃA�E�g�v�b�g */
/* ����(all_sum)��group��(grp_sum) */
%macro eventEvery2weeks(data=, month=, vac_eff=);
  %do i=1 %to %sysfunc(int(((356.25/12) * &month.) / 14 + 1));
    data tmpdata;
	  set &data.;
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
	  label event_sum&i. = "�S�C�x���g���F����&i." event_grp0_sum&i. = "�ΏƌQ�C�x���g���F����&i." event_grp1_sum&i. = "�����Q�C�x���g���F����&i.";
	run;
	/* merged data */
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
  /*
  ods html;
  proc print data=all_sum; run;
  proc print data=grp_sum; run;
  */
  %out_xlsx(data=all_sum, filepath="&dir.\�ʖ�\eventEvery2weeks.xlsx", sheetname="�S�� VacEff: &vac_eff. Month: &month.");
  %out_xlsx(data=grp_sum, filepath="&dir.\�ʖ�\eventEvery2weeks.xlsx", sheetname="�O���[�v�� VacEff: &vac_eff. Month: &month.");
%mend;

/* 
SIM : �V�~�����[�V������
N : �S�ᐔ
ctrl_rate : �ΏƌQ�̔��Ǘ��i%�\���j
vac_eff : vaccine efficacy
month : �ǐՊ���(��)
H0 : �A�����������������ۂ�

e.g.,
SIM = 100, N = 100, ctrl_rate = 1, vac_eff = 0.6, month = 2
==> �V�~�����[�V�����񐔁F100
�S�ᐔ�F100
�ΏƌQ�̔��Ǘ��F1%
vaccine efficacy�F0.6
2�J���Ԃ̒ǐ�
*/
%macro output_Err_Event2weeks(DATA=, SIM=, N=, ctrl_rate=, vac_eff=, month=, alpha=, H0=);
%DataGenerate(DATA=&data., SIM=&SIM., N=&N., ctrl_rate=&ctrl_rate., vac_eff=&vac_eff., month=&month.);
%CalcPvalue(DATA=&data.);
/* ������r */
%calcErr(data=chisq, sim=&SIM., pval=nValue1, alpha=&alpha., H0=&H0.);
%out_xlsx(data=out, filepath="&dir.\�ʖ�\error.xlsx", sheetname="������r�@VacEff: &vac_eff. Month: &month.");
/* ���O�����N���� */
%calcErr(data=logrank, sim=&SIM., pval=ProbChiSq, alpha=&alpha., H0=&H0.);
%out_xlsx(data=out, filepath="&dir.\�ʖ�\error.xlsx", sheetname="���O�����N�@VacEff: &vac_eff. Month: &month.");
/* �|�A�\����A */
%calcErr(data=poisson, sim=&SIM., pval=ProbChiSq, alpha=&alpha., H0=&H0.);
%out_xlsx(data=out, filepath="&dir.\�ʖ�\error.xlsx", sheetname="�|�A�\����A�@VacEff: &vac_eff. Month: &month.");

/* 2�T�Ԗ��̃C�x���g�� */
%eventEvery2weeks(data=&data., month=&month., vac_eff=&vac_eff.);
%mend;

%macro Validate_VacEff(min=, max=, step=, E=, month=);
  %do iter=&min. %to &max. %by &step.;
    %let eff = %sysevalf(&iter./&E.);
    %output_Err_Event2weeks(DATA=scene1, SIM=1000, N=5000, ctrl_rate=0.01, vac_eff=&eff., month=&month., alpha=0.05, H0=FALSE);
  %end;
%mend;

%Validate_VacEff(min=60, max=80, step=5, E=100, month=2);
%Validate_VacEff(min=60, max=80, step=5, E=100, month=2.25);
