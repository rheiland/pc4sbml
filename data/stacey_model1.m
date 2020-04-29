
function driver
%% Settings for solving the differential equations

options = odeset('RelTol',1e-09,...
    'AbsTol',1e-12,...'NonNegative',[1:14],...
    'NormControl','on',...
    'Stats','off',...
    'BDF','off',...
    'MaxOrder', 5);

hours = 24;
tspan = linspace(0,hours*60,100);

% initvalue(1)      = 14;  % MOCCASIN error msg --> 'Identifier' object has no attribute 'lower'

initCond = [14 7.7 1.2 1.9 2.1e-1 0.1 0.88 0.57 0.31 0.061 6.7e-3 0.34 4.9e-2 5.4e-2 8.1 19 0 36]

%% Solve the rate equations

[t, y] = ode15s(@sample_Model_GlycolysisOnly_Core, tspan,initCond,  options);


% ---------------------------------
function dydt = sample_Model_GlycolysisOnly_Core(t,y)

%% Initial Conditions
Glu_in         = y(1);      % Glucose
ATP            = y(2);      % ATP generation
G6P            = y(3);      % Glucose6Phosphate
ADP            = y(4);      % ADP
F6P            = y(5);      % Fructose6Phosphate
FBP            = y(6);      % Fructose16Biphosphate
DHAP           = y(7);      % Dihydroxyacetone phosphate
G3P            = y(8);      % Glyceraldehyde 3-phosphate
NAD            = y(9);      % NAD
NADH           = y(10);     % NADH
thirteenBPG    = y(11);     % 1,3 Biphosphoglycerate
threePG        = y(12);     % 3Phosphoglycerate
twoPG          = y(13);     % 2Phosphoglycerate
PEP            = y(14);     % Phophoenolpyruvate
Pyr            = y(15);     % Pyruvate
Lac_in         = y(16);     % cytoplasmic Lactate
Lac_out        = y(17);     % extracellular Lactate
Glu_out        = y(18);     % extracellular Glucose

%% Glycolysis reactions
%% GLUT1/GLUT3 reaction: Gluout <-> Gluin
Kgluout_glut	= 1.8;
Kgluin_glut     = 10;
Keq_glut        = 1;
f1_glut         = 0.1;
f2_glut         = 0.9;
Keq1_glut       = 1;
Kgluin1_glut	= 10;
Kgluout1_glut	= 9.3;
Vf_glut         = 0.028;

% GLUT_reaction = 1.0
GLUT_reaction = (Vf_glut)*(0.2 + 0.8/(1 - 23.68/50))*((f1_glut*(Glu_out-(Glu_in/Keq_glut)))/(Kgluout_glut*(1 + (Glu_in/Kgluin_glut)) + Glu_out) + ...
               (f2_glut*(Glu_out - (Glu_in/Keq1_glut)))/(Kgluout1_glut*(1 + (Glu_in/Kgluin1_glut)) + Glu_out));
%% HK reaction: Gluin + ATP <-> G6P + ADP
f1_hk       = 0.01;
f2_hk       = 0.99;
Kgluin_hk	= 0.03;
Kg6p_hk     = 0.02;
Katp_hk     = 1.1;
Kadp_hk     = 3.5;
Keq_hk      = 651;
Kgluin2_hk	= 0.3;
Vf_hk       = 0.041;

HK_reaction = (Vf_hk)*(0.2 + 0.8/(1 - 23.68/50))*((((f1_hk/(Kgluin_hk*Katp_hk))*(Glu_in*ATP-(G6P*ADP/Keq_hk)))/...
       (1+ Glu_in/Kgluin_hk + ATP/Katp_hk +((Glu_in*ATP)/(Kgluin_hk*Katp_hk))+ (G6P/Kg6p_hk) + (ADP/Kadp_hk) + ((G6P*ADP)/(Kg6p_hk*Kadp_hk))+((Glu_in*ADP)/(Kgluin_hk*Kadp_hk))+((G6P*ATP)/(Kg6p_hk*Katp_hk)))) ...
      +(((f2_hk/(Kgluin2_hk*Katp_hk))*(Glu_in*ATP-(G6P*ADP/Keq_hk)))/...
      (1 + (Glu_in/Kgluin2_hk) + (ATP/Katp_hk) + ((Glu_in*ATP)/(Kgluin2_hk*Katp_hk)) + ...
       (G6P/Kg6p_hk) + (ADP/Kadp_hk) + ((G6P*ADP)/(Kg6p_hk*Kadp_hk)) + ((Glu_in*ADP)/(Kgluin2_hk*Kadp_hk)) + ((G6P*ATP)/(Kg6p_hk*Katp_hk)))));

%% The HPI reaction : G6P <-> F6P ; E4P, FBP, 6PG

Kg6p_hpi	= 0.4;
Kf6p_hpi	= 0.05;
Kery4p_hpi	= 0.001;
K6pg_hpi	= 0.015;
Kfbp_hpi	= 0.06;
Vf_hpi      = 0.28;
Vr_hpi      = 0.63;

E4P = 0.18;
sixPG = 4.5e-3;
HPI_reaction = ((Vf_hpi)*G6P/Kg6p_hpi - (Vr_hpi)*F6P/Kf6p_hpi)/...
        ( 1 + G6P/Kg6p_hpi + F6P/Kf6p_hpi + E4P/Kery4p_hpi + FBP/Kfbp_hpi + sixPG/K6pg_hpi); 
%% The PFK1 reaction :ATP + F6P <-> FBP + ADP;  F26BP Cit
Kf6p_pfk1	= 1.1;
Kf26bp_pfk1	= 0.00099;
Kiatp_pfk1	= 1.1;
Kcit_pfk1	= 6.7;
alpha_pfk1	= 0.75;
beta_pfk1	= 1.18;
Katp_pfk1	= 0.0292;
L_pfk1      = 6.6;
Kadp_pfk1	= 5;
Keq_pfk1	= 247;
Kfbp_pfk1	= 5;
F26BP       = 0.0034; 
Vf_pfk1     = 0.022;
cCIT = 0.51;

N1_PFK = ((F6P*(1+(F26BP/(alpha_pfk1*Kf26bp_pfk1))))/(Kf6p_pfk1*(1 + (F26BP/Kf26bp_pfk1))));
N2_PFK = (1 + ((F6P*(1+(F26BP/(alpha_pfk1*Kf26bp_pfk1))))/(Kf6p_pfk1*(1+(F26BP/Kf26bp_pfk1)))))^3;
D1_PFK = (L_pfk1*((1+(cCIT/Kcit_pfk1))^4)*((1 + (ATP/Kiatp_pfk1))^4))/((1+(F26BP/Kf26bp_pfk1))^4); 
D2_PFK = ((1+((F6P*(1+(F26BP/(alpha_pfk1*Kf26bp_pfk1))))/(Kf6p_pfk1*(1+(F26BP/Kf26bp_pfk1)))))^4);
PFK1_reaction = (Vf_pfk1)*(0.2 + 0.8/(1 - 23.68/50))*(((ATP/Katp_pfk1)/(1+(ATP/Katp_pfk1)))* ...
                ((1+(beta_pfk1*F26BP/(alpha_pfk1*Kf26bp_pfk1)))/(1+(F26BP/(alpha_pfk1*Kf26bp_pfk1))))* (((N1_PFK*N2_PFK)/(D1_PFK + D2_PFK))...
                - (((ADP*FBP)/(Kadp_pfk1*Kfbp_pfk1*Keq_pfk1))/(ADP/Kadp_pfk1 + FBP/Kfbp_pfk1 + ((ADP*FBP)/(Kadp_pfk1*Kfbp_pfk1)) + 1))));        
%% The ALDO reaction: FBP <-> DHAP + G3P  
Kdhap_aldo	= 0.08;
Kfbp_aldo	= 0.009;
Kg3p_aldo	= 0.16;
Vr_aldo     = 0.063;
Vf_aldo     = 0.08;

ALDO_reaction = (((Vf_aldo)*(FBP/Kfbp_aldo))-((Vr_aldo)*((DHAP*G3P)/(Kdhap_aldo*Kg3p_aldo))))/...
                (1+ FBP/Kfbp_aldo +DHAP/Kdhap_aldo + G3P/Kg3p_aldo +((DHAP*G3P)/(Kdhap_aldo*Kg3p_aldo)));  
 %% The TPI reaction: DHAP <-> G3P
Kdhap_tpi	= 1.6;
Kg3p_tpi	= 0.51;
Vr_tpi      = 28;
Vf_tpi      = 3.4;

TPI_reaction = ((Vf_tpi)*(DHAP/Kdhap_tpi) - (Vr_tpi)*(G3P/Kg3p_tpi))/(1 + DHAP/Kdhap_tpi + G3P/Kg3p_tpi);
%% The GAPDH reaction: NAD + G3P + Pi <-> 13BPG + NADH
Kg3p_gapdh	= 0.19;
Knad_gapdh	= 0.09;
Knadh_gapdh	= 0.01;
Kpi_gapdh	= 11;
Kbpg_gapdh	= 0.022;
cPi         = 7.5; %cPi = 2.5

Vf_gapdh    = 0.331;
Vr_gapdh    = 0.413;

GAPDH_reaction = (((Vf_gapdh)*((NAD*G3P*cPi)/(Knad_gapdh*Kg3p_gapdh*Kpi_gapdh)))-((Vr_gapdh)*((thirteenBPG*NADH)/(Kbpg_gapdh*Knadh_gapdh))))/...
    (1 + NAD/Knad_gapdh +((NAD*G3P)/(Knad_gapdh*Kg3p_gapdh)) + ((NAD*G3P*cPi)/(Knad_gapdh*Kg3p_gapdh*Kpi_gapdh)) + ((thirteenBPG*NADH)/(Kbpg_gapdh*Knadh_gapdh))+(NADH/Knadh_gapdh));
%% The PGK reaction: 13BPG + ADP <-> 3PG + ATP
Ka_pgk      = 0.079;
Kp_pgk      = 0.13;
Kb_pgk      = 0.04;
Kq_pgk      = 0.27;
alpha_pgk	= 1;
beta_pgk	= 1;
Vr_pgk      = 8.7;
Vf_pgk      = 8.7;

PGK_reaction = (((Vf_pgk)*((thirteenBPG*ADP)/(alpha_pgk*Ka_pgk*Kb_pgk)))-((Vr_pgk)*((threePG*ATP)/(beta_pgk*Kp_pgk*Kq_pgk))))/...
    (1 + (thirteenBPG/Ka_pgk) + (ADP/Kb_pgk) + ((thirteenBPG*ADP)/(alpha_pgk*Ka_pgk*Kb_pgk)) + ((threePG*ATP)/(beta_pgk*Kp_pgk*Kq_pgk)) + (threePG/Kp_pgk) + (ATP/Kq_pgk));
%% The PGAM reaction: 3PG <-> 2PG
K3pg_pgam	= 0.19;
K2pg_pgam	= 0.12;
Vr_pgam     = 0.36;
Vf_pgam     = 0.94;

PGAM_reaction = ((Vf_pgam)*threePG/K3pg_pgam - (Vr_pgam)*twoPG/K2pg_pgam)/(1+ threePG/K3pg_pgam+ twoPG/K2pg_pgam);
%% The ENO reaction: 2PG <-> PEP
K2pg_eno        = 0.038;
Kpep_eno        = 0.06;
Vr_eno          = 0.38;
Vf_eno          = 0.34;

ENO_reaction = ((Vf_eno)* twoPG/K2pg_eno - (Vr_eno)*PEP/Kpep_eno)/(1 + twoPG/K2pg_eno + PEP/Kpep_eno);
%% The PYK reaction: PEP + ADP <-> Pyr + ATP
Kadp_pyk        = 0.4;
Katp_pyk        = 0.86;
Keq_pyk         = 1.9517e+05;
Kpep_pyk        = 0.05;
Kpyr_pyk        = 10;
Vf_pyk          = 0.087;

PYK_reaction = Vf_pyk*(((PEP*ADP)/(Kpep_pyk*Kadp_pyk)) - ((Pyr*ATP)/(Kpep_pyk*Kadp_pyk*Keq_pyk)))/((1 + PEP/Kpep_pyk + Pyr/Kpyr_pyk)*(1 + ADP/Kadp_pyk + ATP/Katp_pyk));
%% The LDH reaction: NADH + Pyr <-> Lacin + NAD
Klacin_ldh      = 4.7;
Knad_ldh        = 0.07;
Kpyr_ldh        = 0.3;
Knadh_ldh       = 0.002;
alpha_ldh       = 1;
beta_ldh        = 1;
Vr_ldh          = 0.0740;
Vf_ldh          = 0.0468;

LDH_reaction = (((Vf_ldh)*((NADH*Pyr)/(alpha_ldh*Kpyr_ldh*Knadh_ldh))) - ((Vr_ldh)*((Lac_in*NAD)/(beta_ldh*Klacin_ldh*Knad_ldh))))/...
    (1 + Pyr/Kpyr_ldh + NADH/Knadh_ldh +((NADH*Pyr)/(alpha_ldh*Kpyr_ldh*Knadh_ldh)) + ((Lac_in*NAD)/(beta_ldh*Klacin_ldh*Knad_ldh)) + Lac_in/Klacin_ldh + NAD/Knad_ldh);
%% The MCT reaction: Lac_in <-> Lac_out
Keq_mct1        = 1;
Klacin_mct1     = 8.5;
Klacout_mct1	= 0.5;
Vf_mct1 =	1.45E+00;

MCT_reaction = (Vf_mct1*(Lac_in-(Lac_out/Keq_mct1)))/(Klacin_mct1*(1+(Lac_out/Klacout_mct1))+Lac_in);

%% Rate of Glucose
dydt(1) = GLUT_reaction - HK_reaction;
%% Rate of ATP 
dydt(2) = -HK_reaction - PFK1_reaction + PGK_reaction + PYK_reaction; % + ATP_syn - ATP_util); 
%% Rate of G6P
dydt(3) = HK_reaction - HPI_reaction;
%% Rate of ADP 
dydt(4) = HK_reaction + PFK1_reaction - PGK_reaction - PYK_reaction; % - ATP_syn + ATP_util);  
%% Rate of F6P
dydt(5) = HPI_reaction - PFK1_reaction;
%% Rate of FBP
dydt(6) = PFK1_reaction - ALDO_reaction;
%% Rate of DHAP
dydt(7) = ALDO_reaction - TPI_reaction;
%% Rate of G3P
dydt(8) = ALDO_reaction + TPI_reaction - GAPDH_reaction;
%% Rate of NAD
dydt(9) = -GAPDH_reaction + LDH_reaction;
%% Rate of NADH
dydt(10) = GAPDH_reaction - LDH_reaction;
%% Rate of 13BPG
dydt(11) = GAPDH_reaction - PGK_reaction;
%% Rate of 3PG
dydt(12) = PGK_reaction - PGAM_reaction;
%% Rate of 2PG
dydt(13) = PGAM_reaction - ENO_reaction;
%% Rate of PEP
dydt(14) = ENO_reaction - PYK_reaction;
%% Rate of Pyruvate
dydt(15) = PYK_reaction - LDH_reaction;
%% Rate of Lactate (cytoplasm)
dydt(16) = LDH_reaction - MCT_reaction;
%% Rate of Lactate (extracellular space)
dydt(17) = MCT_reaction;
%% Rate of Glucose (extracellular space)
dydt(18) = -GLUT_reaction;

dydt = dydt'
% return

end

end