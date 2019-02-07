function [Mpp,Mpn,Mnp,Mnn] = bendingReshape(M)

[~,~,C] = size(M);

Mpp = reshape(M(1,1,:),C,1);
Mpn = reshape(M(2,1,:),C,1);
Mnp = reshape(M(1,2,:),C,1);
Mnn = reshape(M(2,2,:),C,1);