function [RNWFE,TNWFE] = AreaChangeBendWFE(E,rho,L,w,wa,ta,beta)

Sa = wa*ta;     %Area of the first section
Sb = Sa*beta;   %Area of the second section

tb = Sb/wa;

Ia = wa*ta^3/12;
Ib = wa*tb.^3/12;

[Ma,Ka] = EB_Beam(rho,Sa,L,E,Ia);
[Mb,Kb] = EB_Beam(rho,Sb,L,E,Ib);

dof = length(Ka)/2;

AA = eye(dof);
BB = eye(dof);
CC = eye(dof);
DD = eye(dof);

%WFE
% 
    nor = 1;
    tol = 1e-4;


    [PhiQ_ApN, PhiQ_AnN, PhiF_ApN, PhiF_AnN]...
      = PolySolve_complex( w,Ka,Ma,L,nor,tol);


    [PhiQ_BpN, ~, PhiF_BpN, ~]...
      = PolySolve_complex( w,Kb,Mb,L,nor,tol);
  
% [PhiQ_ApN, PhiQ_AnN, PhiF_ApN, PhiF_AnN] = WFE_reflection_beamRT_section(Ma,Ka,L,w);
% [PhiQ_BpN, ~, PhiF_BpN, ~] = WFE_reflection_beamRT_section(Mb,Kb,L,w);

%% 
M1 = CC*PhiF_AnN-DD*PhiF_BpN*(BB*PhiQ_BpN)^-1*AA*PhiQ_AnN;
M2 = -CC*PhiF_ApN+DD*PhiF_BpN*(BB*PhiQ_BpN)^-1*AA*PhiQ_ApN;

N1 = DD*PhiF_BpN-CC*PhiF_AnN*(AA*PhiQ_AnN)^-1*BB*PhiQ_BpN;
N2 = CC*PhiF_ApN-CC*PhiF_AnN*(AA*PhiQ_AnN)^-1*AA*PhiQ_ApN;

RNWFE = M1\M2;
TNWFE = N1\N2;


end

