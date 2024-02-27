function fac = normint(c)

% calculates normalization integrals of modes

% define coefficients for +/- integrals
Cpa = 1./sqrt(2)./c.kperpi.*(-i.*c.beta.*c.CmE + c.k.*c.mu.*c.CmH);
Cma = 1./sqrt(2)./c.kperpi.*(i.*c.beta.*c.CmE + c.k.*c.mu.*c.CmH);
Cp = -1./sqrt(2)./c.kperpi.*(i.*c.beta.*c.CmE + c.k.*c.mu.*c.CmH);
Cm = 1./sqrt(2)./c.kperpi.*(i.*c.beta.*c.CmE - c.k.*c.mu.*c.CmH);

% analytical integral components
% check minus signs
analEzEz = c.CmE*c.CmE*bessint(c.orders, c.kperpi, c.a);
analEpEm = Cpa.*Cm.*bessint(c.orders-1, c.kperpi, c.a);
analEmEp = Cma.*Cp.*bessint(c.orders+1, c.kperpi, c.a);

% total
tot = analEzEz + analEpEm + analEmEp;
fac = sqrt(tot);

function val = bessint(m, k, a)
% evaluates Bessel integral
% check minus sign
val = pi.*a.^2.*(besselj(m, k.*a).^2 - besselj(m-1, k.*a).*besselj(m+1, k.*a));
