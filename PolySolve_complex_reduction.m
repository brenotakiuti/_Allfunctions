function [ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,llRpa,llRna,sa,kpa,kna,Phiq,pureN] = PolySolve_complex_reduction( w,Ka,Ma,La,nor,tol,l_rm,lim,lim2)
%% Using fan2016 and fan2018* as reference
%02/01/2019
%
% *Y. Fan, C.W. Zhou, J.P. Laine, M. Ichchou, L. Li,Model reduction schemes 
% for the wave and finite element method using the free modes of a unit 
% cell, Computers & Structures,Volume 197,2018.
%
% l_rm: number of modes retained
%
% lim: used in pi/lim to calculate the tolerance to distinguish pure real
% and pure imaginary. bigger lim = tighter tolerance; smaller lim = looser
% tolerance
%
% lim2: used to distinguish between useful results and numerical modes.
% smaller is tighter
%
% tol: used to distinguish zeroes from very small values. Smaller is
% tighter
%
% IMPORTANT: when using this function, the PSI (left eigenvector) must be
% calculated for each frequency and used to calculate the scattering
% matrix. In this case using the lowest frequency PSI for all frequencies
% does not work.
l_rm = 82;

if nargin<8
    lim = 5;
    lim2 = 2000;
end

    [dofa,~] = size(Ka);
    dofa = dofa/2;

% PHIq and PHIf for the first section (a)
    % Define the dinamic matrix
    D = (Ka-w^2*Ma);
%     D = (Ka-w(q)^2*Ma);
    
    % Saparate it in quadrants
    DLL = (D(1:dofa,1:dofa));
    DLR = (D(1:dofa,dofa+1:dofa*2));
    DRL = (D(dofa+1:dofa*2,1:dofa));
    DRR = (D(dofa+1:dofa*2,dofa+1:dofa*2));
    
    Z=zeros(max(size(DRR)));
    I=eye(max(size(DRR)));
%     
%     A2 = DLR;
%     A1 = DLL+DRR;
%     A0 = DRL;

%     DD=[Z,I;-pinv(DLR)*DRL,-pinv(DLR)*(DLL+DRR)];
    
%     DD=[Z,I;-inv(DLR)*DRL,-inv(DLR)*(DLL+DRR)];   %changed this because
%     of the hybrid version. If it starts giving strange results. Restore.

    % Zhong
%     DD = [ -pinv(DLR)*DLL pinv(DLR); -DRL+DRR*pinv(DLR)*DLL -DRR*pinv(DLR)];
    %*Dont forget to change the Phif later
    
    % FAN and ICHCHOU 2016
    [~,N2] = size(DRR);
    N2 = N2^2;
    sig = norm(DRR,2)/N2;
    
    AA = [ Z sig*I; -DRL -DRR];
    BB = [ sig*I Z; DLL DLR];
    
    [VR,llR] = eig(AA,BB);
%     VL_nmlzr = inv(VR*VL);
%     VLn = VL*VL_nmlzr;
%     [VR,llR] = eigs(DD);
%     [VR,llR]=polyeig(A0,A1,A2);
    
    sa=1i.*log((diag(llR)));   %%llR=exp(-i*s*L) I find s
    sa=sa./La; %real positive value are for positive propagating
    % ee=10^(-4); %error
%     ee = tol;
    
    %% Normalize PhiQ!!!!
    
    [r,~] = size(VR);
    

    Phiq = zeros(size(VR));
    switch nor
        case {1,2}           
            for ii = 1:r
                Phiq(ii,:) = VR(ii,:)./(VR(nor,:));
            end
        case 3
            for ii = 1:r
                Phiq(ii,:) = VR(ii,:)./max(VR);
            end
        otherwise
            Phiq = VR;
    end
    llRd = diag(llR);
    
   
    %% First simply organize the vectors
    % Here we simple sort the vectors by the absolute value of k
    
    [~,index] = sort(abs(sa));
    
    % Join the sorted vector with its indexing vector
    s1 = [sa(index) index];
    
    % Reproduce the sorting order into Phiq
    Phiq1 = Phiq(:,index);

    %% Mode reduction by fan2018
    %   There are 5 types of model reduction cited by Fan2018.
    %   Here we try all of then, and comparison will tell which is 
    % better for our application. All of then differ only in the way
    % the residual flexibility is used. The low order H (H_low) is
    % calculated in the same way for all cases.
    
%     H_low = 
    
    % Model reduction 1: Free(0th mode)
        
end