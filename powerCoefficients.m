function [PR,PT] = powerCoefficients(TRT,Pa,Pb)

[na,~] = size(Pa);
na = na/2;
[nb,~] = size(Pb);
nb = nb/2;

for ii = 1:na
    for jj = 1:na
        PR(jj,ii)=abs(TRT(jj,ii))^2*(Pa(na+jj,na+jj)/Pa(ii,ii));
    end
end

for ii = 1:na
    for jj = 1:nb
        PT(jj,ii)=abs(TRT(na+jj,ii))^2*(Pb(jj,jj)/Pa(ii,ii));
    end
end

%     PrPP2(q) = abs(RPAA(1,1,q))^2*(Pa2(nmodes_a+1,nmodes_a+1)/Pa2(1,1));
%     PrPL2(q) = abs(RPAA(2,1,q))^2*(Pa2(nmodes_a+2,nmodes_a+2)/Pa2(1,1));
%     PrPN2(q) = abs(RPAA(3,1,q))^2*(Pa2(nmodes_a+3,nmodes_a+3)/Pa2(1,1));
% %     PrPC12(q) = abs(RPAA(4,1,q))^2*(Pa2(nmodes_a+4,nmodes_a+4)/Pa2(1,1));
% %     PrPC22(q) = abs(RPAA(5,1,q))^2*(Pa2(nmodes_a+5,nmodes_a+5)/Pa2(1,1));
%     PtPP2(q) = abs(TPCA(1,1,q))^2*(Pc2(1,1)/Pa2(1,1));
%     PtPL2(q) = abs(TPCA(2,1,q))^2*(Pc2(2,2)/Pa2(1,1));
%     PtPN2(q) = abs(TPCA(3,1,q))^2*(Pc2(3,3)/Pa2(1,1));
% %     PtPC12(q) = abs(TPCA(4,1,q))^2*(Pc2(4,4)/Pa2(1,1));
% %     PtPC22(q) = abs(TPCA(5,1,q))^2*(Pc2(5,5)/Pa2(1,1));
%    
%     PrLP2(q) = abs(RPAA(1,2,q))^2*(Pa2(nmodes_a+1,nmodes_a+1)/Pa2(2,2));
%     PrLL2(q) = abs(RPAA(2,2,q))^2*(Pa2(nmodes_a+2,nmodes_a+2)/Pa2(2,2));
%     PrLN2(q) = abs(RPAA(3,2,q))^2*(Pa2(nmodes_a+3,nmodes_a+3)/Pa2(2,2));
% %     PrLC12(q) = abs(RPAA(4,2,q))^2*(Pa2(nmodes_a+4,nmodes_a+4)/Pa2(2,2));
% %     PrLC22(q) = abs(RPAA(5,2,q))^2*(Pa2(nmodes_a+5,nmodes_a+5)/Pa2(2,2));
%     PtLP2(q) = abs(TPCA(1,2,q))^2*(Pc2(1,1)/Pa2(2,2));
%     PtLL2(q) = abs(TPCA(2,2,q))^2*(Pc2(2,2)/Pa2(2,2));
%     PtLN2(q) = abs(TPCA(3,2,q))^2*(Pc2(3,3)/Pa2(2,2));
% %     PtLC12(q) = abs(TPCA(4,2,q))^2*(Pc2(4,4)/Pa2(2,2));
% %     PtLC22(q) = abs(TPCA(5,2,q))^2*(Pc2(5,5)/Pa2(2,2));