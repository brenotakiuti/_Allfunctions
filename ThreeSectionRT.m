function [RAA,TCA] = ThreeSectionRT(Raa_1,Rbb_2,Rbb_1,Tba_1,Tcb_2,Tab_1,Tb,Tbc_2,Rcc_2)

% Reflection coefficients and phase change:
% Beam in logitudinal motion with 3 section discontinuity
% Breno Ebinuma Takiuti
% 28/08/2016

[dofb,~]=size(Rbb_1);

RAA = Raa_1 + Tab_1*pinv(eye(dofb)-(Tb)*Rbb_2*Tb*Rbb_1)*Tb*Rbb_2*Tb*Tba_1;
TCA = Tcb_2*pinv(eye(dofb)-Tb*Rbb_1*(Tb)*Rbb_2)*(Tb)*Tba_1;

% RCC = Rcc_2 + Tcb_2*pinv(eye(dofb)-Tb*Rbb_1*pinv(Tb)*Rbb_2)*pinv(Tb)*Rbb_1*pinv(Tb)*Tbc_2;
% TAC = Tab_1*pinv(eye(dofb)-pinv(Tb)*Rbb_2*Tb*Rbb_1)*Tb*Tbc_2;