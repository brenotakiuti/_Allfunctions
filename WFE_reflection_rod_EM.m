% clear
% clc 
% 
% %% Material's Constants
% % Material: Aluminium
% 
% E = 70e9;               % Young's Modulus (Pa)
% rho = 2.7e3;            % density (kg/m^3)
% 
% %% Geometric constants
% 
% A = (50e-3)^2;          % Cross-sectional area (m^2)
% L = 0.01;                % length of the element (m)
% 
% %% M & K
% % Defined by the finite element analysis
% 
% K = E*A/L*[1 -1; -1 1];
% M = rho*A*L*[2/6 1/6; 1/6 2/6];
% w = 100;

%% 
function [PhiQ_pL, PhiQ_nL, PhiF_pL, PhiF_nL,PhiQ_pR, PhiQ_nR, PhiF_pR, PhiF_nR, kpp, kpn] = WFE_reflection_rod_EM(M,K,L,w)

%% M & K
% Defined by the finite element analysis

dof = length(M)/2;  % total degrees of freedom in one side of the element

%% Eigenvalue problems
%% Polynomial eigenvalue problem
% known frequency vector (omega), find k

% Define the dinamic matrix
D = (K-w^2*M);

% Saparate it in quadrants
DLL = (D(1:dof,1:dof));
DLR = (D(1:dof,dof+1:dof*2));
DRL = (D(dof+1:dof*2,1:dof));
DRR = (D(dof+1:dof*2,dof+1:dof*2));

% Define the polinomial eigenvalue problem
%   (lambda^2*DLR + lambda*(DLL+DRR) + DRL)qL = 0
%   (lambda^2*A2  + lambda*A1        + A0 )qL = 0
A2 = DLR;
A1 = DLL+DRR;
A0 = DRL;

%% Other formulation, Renno (2013)
% DD = [ -inv(DLR)*DLL inv(DLR); -DRL+DRR*inv(DLR)*DLL -DRR*inv(DLR)];    % Renno(2013)

%% Solve the eigenvalue problem
[Phi,lam] = polyeig(A0,A1,A2);

ke = roundn((-log(lam)/L),-5);  % caution with truncation errors
% ke = (-log(lam)/L);

%% Normalize PhiQ!!!!

[r,c] = size(Phi);

Phiq = Phi;
for ii = 1:r
    Phiq(ii,:) = Phi(ii,:)./Phi(1,:);
end

%% Order the eigen solutions

pp = 0; pn = 0; 

for ii = 1:length(ke)
    Phif(:,ii) = D*[Phiq(:,ii);lam(ii)*Phiq(:,ii)];
    if isreal(lam(ii))==0 %propagating
        if imag(ke(ii))>0
            
            pp = pp+1;
            kpp(pp) = ke(ii);  %positive propagating
            Phiq_pp(:,pp) = Phiq(:,ii);
            Phif_pp(:,pp) = Phif(1:dof,ii);
            Phiq_ppR(:,pp) = lam(ii)*Phiq(:,ii);
            Phif_ppR(:,pp) = Phif(dof+1:2*dof,ii);
            
        elseif imag(ke(ii))<0
            
            pn = pn+1;
            kpn(pn) = ke(ii);  %negative propagating
            Phiq_pn(:,pn) = Phiq(:,ii);
            Phif_pn(:,pn) = Phif(1:dof,ii);
            Phiq_pnR(:,pn) = lam(ii)*Phiq(:,ii);
            Phif_pnR(:,pn) = Phif(dof+1:2*dof,ii);
            
        end
    end
end

PhiQ_pL = [Phiq_pp];

PhiQ_nL = [Phiq_pn];

PhiF_pL = [Phif_pp];

PhiF_nL = [Phif_pn];

PhiQ_pR = [Phiq_ppR];

PhiQ_nR = [Phiq_pnR];

PhiF_pR = [Phif_ppR];

PhiF_nR = [Phif_pnR];

% end

