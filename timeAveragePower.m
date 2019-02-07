function P = timeAveragePower(AMP,Pa)
[na,~] = size(Pa);
P = zeros(na,1);
for ii = 1:na
    
    AMP_I = zeros(size(AMP));
    AMP_I(ii) = AMP(ii);

    P(ii) = (AMP_I'*Pa*AMP_I);
    
end