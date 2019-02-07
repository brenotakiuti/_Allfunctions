function [PhiQ_pT,PhiQ_nT,PhiF_pT,PhiF_nT] = WM_reflection_rod(rho,A,E,w)

k = abs(sqrt(rho/E)*w);

PhiQ_pT = 1;

PhiQ_nT = 1;

PhiF_pT = -1i*E*A*k;

PhiF_nT = 1i*E*A*k;

