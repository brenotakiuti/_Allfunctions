function [ni,tau,PWi,PWr,PWt] = powerRatio(w, PhiQa_p, PhiQa_n, PhiQc_p, PhiQc_n, PhiFa_p, PhiFa_n, PhiFc_p, PhiFc_n, AMP_i, AMP_r, AMP_t)
% Calculate the power coefficients by dividing incident power into the
% reflected and transmitted powers
% 

qi = [PhiQa_p PhiQa_n]*AMP_i;
qr = [PhiQa_p PhiQa_n]*AMP_r;
qt = [PhiQc_p PhiQc_n]*AMP_t;
fi = [PhiFa_p PhiFa_n]*AMP_i;
fr = [PhiFa_p PhiFa_n]*AMP_r;
ft = [PhiFc_p PhiFc_n]*AMP_t;

% Incident Power
PWi = -w/2*imag(fi'*qi);
% Reflected Power
PWr = -w/2*imag(fr'*qr);
% Transmitted Power
PWt = -w/2*imag(ft'*qt);

ni = PWr/PWi;
tau = PWt/PWi;