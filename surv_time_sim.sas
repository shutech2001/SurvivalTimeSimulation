%let dir = C:\Users\Shimizu\���������J���@�l �������ۈ�Ì����Z���^�[\biostat share - General\Mpox_RCT;

ods html close;

/*
�V�~�����[�V�����f�[�^�����p�̃}�N��

DATA : �V�~�����[�V�����p�ɍ쐬����f�[�^��
SIM : �V�~�����[�V������
N : �ᐔ�i���Q���킹�����́D�������e�Q�ɕ�������D�j
CONTROL_PROP : �ΏƌQ�̔��Ǌm��
VE : vaccine efficacy
MONTH : �ώ@����
*/
%macro DataGenerate(DATA=, SIM=, N=, CONTROL_PROP=, VE=, MONTH=);
data &DATA.;
  call streaminit(42);
  do SIM = 1 to &SIM.;
    do ID = 1 to &N.;
      CONTROL_PROP = &CONTROL_PROP.; VE = &VE.;
	  /* �ΏƌQ */
      if ID <= &N./2 then do;
	    GRP = 0; P = CONTROL_PROP; end;
	  else do;
	    GRP = 1; P = (1 - VE) / 100; end;
      /* �ݐϊm�� = p ��� 1/lambda���v�Z */
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
P�l���v�Z���邽�߂̃}�N��
���X�N���E���O�����N����E�|�A�\����A��3�̉�͂��s��

DATA : ��͑ΏۂƂȂ�f�[�^�Z�b�g��
*/
%macro CalcPvalue(DATA=);
  /* ���X�N�� */
  proc freq data = &DATA.;
    by sim;
    table grp*event / riskdiff(equal column=1 var=sample) chisq nocol nopercent;
    ods output PdiffTest=chisq(where=(Name1="P2_RDIF1"));
  run;
  /* ���O�����N */
  proc lifetest data = &DATA.;
    by sim;
    time time*censor(1);
    strata / group = grp;
    ods output HomTests=logrank(where=(Test="���O�����N"));
  run;
  /* �|�A�\����A */
  proc genmod data=&DATA.;
    by sim;
    class grp;
    model event(ref="0") = grp / dist = poisson link = log offset = TIME;
    ods output ParameterEstimates=poisson(where=(Parameter="GRP" and Level1 = "0"));
  run;
%mend;

/*
alpha error, beta error���v�Z���邽�߂̃}�N��
�o�̓f�[�^�̓G���[�̒l��1�����i�[���ꂽ�f�[�^�ƂȂ�

DATA : ��͌��ʂ��i�[���ꂽ�f�[�^
OUT : �G���[�����i�[�����o�̓f�[�^�̖��O
SIM : �V�~�����[�V�����񐔁i�G���[���Z�o���邽�߂ɕK�v�j
PVAL : ��͌��ʃf�[�^��P�l�̗�̖��O
ALPHA : �L�Ӑ���
H0 : H0�����������ۂ��iif H0 is TRUE then calculate alpha error, otherwise calculate beta error (power)�j
*/
%macro calcErr(DATA=, OUT=, SIM=, PVAL=, ALPHA=, H0=);
  data &DATA.;
    set &DATA.;
	if &PVAL. < &ALPHA. then rej_flg=1; else rej_flg=0; run;
  /* �G���[�t���O���Ă� */
  data flg;
    set &DATA.;
	retain rej 0; rej+rej_flg; *reject�𐔂��グ�� (for alpha error);
    retain acc 0; acc+abs(1-rej_flg); *accept�𐔂��グ�� (for beta error);
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
DATA : �f�[�^
FILEPATH : �t�@�C���̃p�X
SHEETNAME : �V�[�g�̖��O
*/
%macro out_xlsx(DATA=, FILEPATH=, SHEETNAME=);
  proc export data = &data. outfile = &filepath. dbms = xlsx replace;
    sheet = &sheetname.;
    label;
  run;
%mend;

/*
�o�͌`���𐮂����e�[�u�����o�͂���}�N��

DATA : �V�~�����[�V�����p�f�[�^�̖��O
SIM : �V�~�����[�V������
N : �ᐔ�i���Q���킹�����́j
CONTROL_PROP : �ΏƌQ�̔��Ǌm��
VE : vaccine efficacy
MONTH : �ώ@���ԁi���j
H0 : H0�����������ۂ�
ALPHA : �L�Ӑ���
SET : ���łɂ���f�[�^�Z�b�g�ɏc�������邩�ۂ��i2�ڂ̃V�i���I�ȍ~�͂��Ƃ��Ƃ���e�[�u���ɏc�������Ă����j
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
  /* �c���� */
  %if &SET. = TRUE %then %do;
    data outTab; set outTab partTab; run;
  %end;
  %if &SET. = FALSE %then %do;
    data outTab; set partTab; run;
  %end;
%mend;
/*
2�T�Ԗ��ɃC�x���g�������������Ă��邩������悤�ȃA�E�g�v�b�g
����(all_sum)��group��(grp_sum) 

DATA : �V�~�����[�V�����p�̃f�[�^�Z�b�g��
MONTH : �ώ@����
CONTROL_PROP : �ΏƌQ�̔��Ǌm��(�V�[�g���ɋL�ڂ��邽��)
VE : vaccine efficacy�i�V�[�g���ɋL�ڂ��邽�߁j
*/
%macro eventEvery2weeks(DATA=, MONTH=, CONTROL_PROP=, VE=);
  %do i=1 %to %sysfunc(int(((356.25/12) * &MONTH.) / 14 + 1));
    data tmpdata;
	  set &DATA.;
      /* 2�T�Ԃ��؂�o�� */
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
    /* 2���Ԗڈȍ~�͉����� */
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
  %out_xlsx(DATA=all_sum, FILEPATH="&dir.\�ʖ�\eventEvery2weeks.xlsx", SHEETNAME="�S�� CONTROL_PROP = &CONTROL_PROP. VE = &VE. MONTH = &MONTH.");
  %out_xlsx(DATA=grp_sum, FILEPATH="&dir.\�ʖ�\eventEvery2weeks.xlsx", SHEETNAME="�O���[�v�� CONTROL_PROP = &CONTROL_PROP. VE = &VE. MONTH = &MONTH.");
%mend;

/*
CONTROL_PROP, VE�̒l�𓮓I�ɕς��ăV�~�����[�V���������s���邽�߂̃}�N��

%do�}�N���̎d�l��C�X�e�b�v����int�^�łȂ��Ƃ����Ȃ����ߎg�����v����

e.g., CONTROL_PROP��0.005 ~ 0.010�܂�0.001���ς��ăV�~�����[�V�������s�������Ƃ�
  ==> CTRL_MIN = 5, CTRL_MAX = 10, CTRL_STEP = 1, CTRL_E = 1000
         i.e., �S�Đ����l�ɂȂ�悤��CTRL_E�{����iVE�̏ꍇ�����l�j

<< �V�~�����[�V�����񐔓���ύX����ꍇ�́C�}�N������SIM����ύX >>

CTRL_MIN : CONTROL_PROP�̍ŏ��l / CTRL_E
CTRL_MAX : CONTROL_PROP�̍ő�l / CTRL_E
CTRL_STEP : CONTROL_PROP�̃X�e�b�v�� / CTRL_E
CTRL_E : 
VE_MIN : VE�̍ŏ��l / VE_E
VE_MAX : VE�̍ő�l / VE_E
VE_STEP : VE�̃X�e�b�v�� / VE_E
VE_E : 
MONTH : �ώ@���ԁi�V�[�g���Ƃ��Ă���j
*/
%macro Validate_CTRL_VE(CTRL_MIN=, CTRL_MAX=, CTRL_STEP=, CTRL_E=, VE_MIN=, VE_MAX=, VE_STEP=, VE_E=, MONTH=);
  %do p_ = &CTRL_MIN. %to &CTRL_MAX. %by &CTRL_STEP.;
    %let p = %sysevalf(&p_./&CTRL_E.);
    %do ve_ = &VE_MIN. %to &VE_MAX. %by &VE_STEP.;
      %let ve = %sysevalf(&ve_./&VE_E.);
      %if &p_. = &CTRL_MIN. and &ve_. = &VE_MIN. %then %createErrorTab(DATA=data, SIM=1000, N=5000, CONTROL_PROP=&p., VE=&ve., MONTH=&MONTH., H0=FALSE, ALPHA=.05, SET=FALSE);
      %else %createErrorTab(DATA=data, SIM=1000, N=5000, CONTROL_PROP=&p., VE=&ve., MONTH=&MONTH., H0=FALSE, ALPHA=.05, SET=TRUE);
      /* 2�T���̃C�x���g�����o���Ƃ��͉��̃R�����g�A�E�g�����s�i���s���Ԃ�����j */
      *%eventEvery2weeks(DATA=data, MONTH=&MONTH., CONTROL_PROP=&p., VE=&ve.);
    %end;
  %end;
  %out_xlsx(data=outTab, filepath="&dir.\�ʖ�\error.xlsx", sheetname="MONTH=&MONTH.");
  /* outTab���폜���Ĉӎ������Ɋ��ɂ���f�[�^�ɏc�������邱�Ƃ�h�� */
  proc delete data=outTab; run;
%mend;

/* CONTROL_PROP = 0.01, VE = 0.5~0.8�̗� */
%Validate_CTRL_VE(CTRL_MIN=1, CTRL_MAX=1, CTRL_STEP=1, CTRL_E=100, VE_MIN=5, VE_MAX=8, VE_STEP=1, VE_E=10, MONTH=2);
