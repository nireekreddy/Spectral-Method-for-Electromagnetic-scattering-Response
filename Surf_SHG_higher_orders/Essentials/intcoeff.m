function Cm = intcoeff(order, Am, Bm, kperp, kperpi, a)

% Interior Bessel coefficients from exterior coefficients

nka = kperp.*a;
nika = kperpi.*a;

Cm = (Am.*besselj(order, nka) + Bm.*besselh1(order, nka))./besselj(order, nika);
