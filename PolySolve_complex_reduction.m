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
    
   
    
    %% S1: Separate only the pure real and pure imaginary wavenumbers
%     oc = unique_tol(sa,tol);
%     orderV = [];
    count = 1;
%     tol2 = 1e-4;
    tol2 = tol;
    
    % Generate a matrix of a combination of wavenumbers and their original
    % position on the wavenumbers vector: [k 1:1:end]_MODESx2
    s1 = [sa (1:length(sa))'];
    % Obs.: this is a trick to integrate the vector indexes into the
    % sorting process, so that any change in the vector's position is
    % "recorded" along the process. With this the second column of the
    % matrix can be used to "repeat the sorting" of the wavenumbers to the
    % wavemodes, without having to sort both along the process.
    
    s2 = s1;
    for ii = 1:length(sa)

%         if (abs(real(sa(ii)*La))<pi/lim && abs(imag(sa(ii)*La))<pi/lim)
%             % Separate the numbers lower than pi over something
            if (abs(real(sa(ii)))>tol2 && abs(imag(sa(ii)))<tol2) || (abs(real(sa(ii)))<tol2 && abs(imag(sa(ii)))>tol2)
                % Separate pure real and pure imag and put then in the
                % beginning of the vector
                s2(count,:) = s1(ii,:);
                s2(ii,:) = s1(count,:);
                count = count+1;
            end
%         end
    end
    nSort = count-1;

    % Sort the pure numbers:
    % First real then imaginary ('descend')
    warning('off','all')
%     s2
%     nSort
    [Pos, Neg] = pureRISort(s2(1:nSort,:),'descend');
    warning('on','all')
    
    % Idea! Should we use the l_rm based on the size of Pos?    size(s2)
    size(Pos);  
    
%     size(Pos)
    % Remove the sorted number, retain only the complex (including errors)
    s5 = s2;
    s5(1:nSort,:) = [];

    
    %% Take only COMPLEX modes
    
    orderV = [];
    s8 = [];
    count = 1;
%     tol2 = 1e-4;
    tol2 = tol;
    for ii = 1:length(s5)
        if abs(real(s5(ii)))<lim2 && abs(real(s5(ii)))>tol2
            % this value of 400 is like pi/6: is an arbitrary value that
            % functions a separation between numerical errors and the valid results 
            s8(count,:) = s5(ii,:);
            orderV(count)=ii;
            count = count+1;
        end
    end
%     s8;
    % Remove the complex small numbers, retain the rest
    s5(orderV,:) = [];

    s8 = s5;

    % Sort the complex numbers: first the real positive, if equal real,
    % imaginary negative first
    warning('off','all')
    [Cpos, Cneg] = complexSort2(s8,'ascend',1e-2);
%     [Cpos, Cneg] = complexSort2(s5,'ascend',1e-2);      %Just use everything!
%     [Cpos, Cneg] = complexSort2(s8,'descend',1e-2);
    warning('on','all')
    size(Cpos)
    size(Cneg)
    %% Join back all the numbers:
    % First pure real, second the pure imaginary, third the complex and
    % last the errors
%     nRest = length(s5(:,1));
    
    % Sort the errors in the same way as the complex numbers, without
    % mixing with the other 3 types of numbers
%     [rp,rn] = complexSort(s5,'ascend');
    
    % Rebuild the wavenumbers vectors: separated in positive and negative
    % going waves
%     s6p = [Pos; Spos; Cpos; rp];
%     s6n = [Neg; Sneg; Cneg; rn];
%     s6p = [Pos];
%     s6n = [Neg];
    pureN = length(Pos);
    s6p = [Pos; Cpos];
    s6n = [Neg; Cneg];   
 
%     s6p = [Pos; Cpos; rp];
%     s6n = [Neg; Cneg; rn];

    % Separate the ordering vector
    index_p = real(s6p(:,2));
    index_n = real(s6n(:,2));
    
    [ppa,~] = size(s6p);
    [nn,~] = size(s6n);
    
    % Separate the wavenumbers vectors
    kpa = s6p(1:ppa,1);
    kna = s6n(1:nn,1);
    
    % Separate the lambda
    llRpa = llRd(index_p);
    llRna = llRd(index_n);
    
    % Calculate the "unsorted" wavemodes and left eigenvectors
%     Phif = D*Phiq;
%     Psiq = pinv(Phiq);  %Another way is to calculate VR*inv(VR*VL) usinge VL from eig or VL'=eig(AA')
%     Psif = pinv(Phif);
    
    Phiq2 = Phiq(:,[index_p index_n]);
%     Phif2 = Phiq(:,[index_p index_n]);
    
    Phiql1 = Phiq2(1:dofa,:);
%     Phifa1 = Phif(1:dofa,:);
    Phiqr1 = Phiq2(dofa+1:2*dofa,:);     %Zhong
    Phifl1 = DLL*Phiql1 + DLR*Phiqr1;
%     Phifr1 = DRL*Phiql1 + DRR*Phiqr1;

    PhiQa_p = Phiql1(:,1:length(index_p));
    PhiQa_n = Phiql1(:,length(index_p)+1:end);
    PhiFa_p = Phifl1(:,1:length(index_p));
    PhiFa_n = Phifl1(:,length(index_p)+1:end);
%%   MITROU  
%     PHI = [PhiQa_p PhiQa_n;PhiFa_p PhiFa_n];
%     PSI = pinv(PHI);
%     
%     Psiq = PSI(:,1:dofa);
%     Psif = PSI(:,dofa+1:end);
%     
%     PsiQa_p = Psiq(1:length(index_p),:);
%     PsiQa_n = Psiq(length(index_p)+1:end,:);
%     PsiFa_p = Psif(1:length(index_p),:);
%     PsiFa_n = Psif(length(index_p)+1:end,:);
%%  Another option (FAN/ICHCHOU)

    PHIp = [PhiQa_p ;PhiFa_p ];
    PSIp = pinv(PHIp);
    PHIn = [PhiQa_n ;PhiFa_n ];
    PSIn = pinv(PHIn);
    
%     PsiQa_p = PSIp(:,1:dofa);
%     PsiFa_p = PSIp(:,dofa+1:end);
%     PsiQa_n = PSIn(:,1:dofa);
%     PsiFa_n = PSIn(:,dofa+1:end);
    
    PsiFa_p = PSIp(:,1:dofa);
    PsiQa_p = PSIp(:,dofa+1:end);
    PsiFa_n = PSIn(:,1:dofa);
    PsiQa_n = PSIn(:,dofa+1:end);

    %% Mode reduction by fan2018
    %   There are 5 types of model reduction cited by Fan2018.
    %   Here we try all of then, and comparison will tell which is 
    % better for our application. All of then differ only in the way
    % the residual flexibility is used. The low order H (H_low) is
    % calculated in the same way for all cases.
    
%     H_low = 
    
    % Model reduction 1: Free(0th mode)
        
end

