function [RNWFE,TNWFE] = AreaChangeLongWFE(E,rho,L,w,Sa,beta)

% Eingenproblem solution of a Rod (2 nodes, 1dof per node)
% Area change problem
% Breno Ebinuma Takiuti
% 14/10/2016

%% Solution parameters

% Area
Sb = Sa*beta;

%% M & K
% Defined by the finite element analysis

Ka = E*Sa/L*[1 -1; -1 1];
Ma = rho*Sa*L*[2/6 1/6; 1/6 2/6];
    
Kb = E*Sb./L*[1 -1; -1 1];
Mb = rho*Sb.*L*[2/6 1/6; 1/6 2/6];

dof = length(Ka)/2;

%% WFE

[PhiQ_ApN, PhiQ_AnN, PhiF_ApN, PhiF_AnN] = WFE_reflection_rod(Ma,Ka,L,w);
[PhiQ_BpN, PhiQ_BnN, PhiF_BpN, PhiF_BnN] = WFE_reflection_rod(Mb,Kb,L,w);

A1 = eye(dof);
A2 = eye(dof);
B1 = eye(dof);
B2 = eye(dof);

Cn = -B1*PhiF_AnN + B2*PhiF_BpN*(A2*PhiQ_BpN)^-1*A1*PhiQ_AnN;
Dn = B1*PhiF_ApN - B2*PhiF_BpN*(A2*PhiQ_BpN)^-1*A1*PhiQ_ApN;
En = -B2*PhiF_BpN + B1*PhiF_AnN*(A1*PhiQ_AnN)^-1*A2*PhiQ_BpN;
Fn = -B1*PhiF_ApN + B1*PhiF_AnN*(A1*PhiQ_AnN)^-1*A1*PhiQ_ApN;

RNWFE = roundn((inv(Cn)*Dn),-4);
TNWFE = roundn((inv(En)*Fn),-4);

% T = [A1*PhiQ_AnN -B1*PhiQ_BpN; A2*PhiF_AnN -B2*PhiF_BpN]^-1*[-A1*PhiQ_ApN B1*PhiQ_BnN; -A2*PhiF_ApN B2*PhiF_BnN];


