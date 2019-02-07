function PhiO = unwrap2(PhiI,tol)

for ii=1:length(PhiI)-1
%     abs(PhiI(ii))-abs(PhiI(ii+1));
        abs(abs(PhiI(ii))-abs(PhiI(ii+1)))
    while abs(abs(PhiI(ii))-abs(PhiI(ii+1)))>=tol
        PhiI(ii) > PhiI(ii+1);

        if PhiI(ii) < PhiI(ii+1)
            PhiI(ii+1:end) = PhiI(ii+1:end)-pi;
        else 
            PhiI(ii+1:end) = PhiI(ii+1:end)+pi;
        end
        abs(abs(PhiI(ii))-abs(PhiI(ii+1)))
    end
end

PhiO = PhiI;