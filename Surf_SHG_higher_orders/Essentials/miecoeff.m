function [Mee, Meh, Mhe, Mhh] = miecoeff(c)

Jm = besselj(c.orders, c.kperp.*c.a); dJm = besseljd(c.orders, c.kperp.*c.a);
Hm = besselh(c.orders, c.kperp.*c.a); dHm = besselh1d(c.orders, c.kperp.*c.a);
Jmi = besselj(c.orders, c.kperpi.*c.a); dJmi = besseljd(c.orders, c.kperpi.*c.a);

e1 = -c.orders.*c.beta./c.a.*Jm.*(1./c.kperpi.^2 - 1./c.kperp.^2);
e2 = i.*c.k.*(-c.mui./c.kperpi.*dJmi./Jmi.*Jm + c.mu./c.kperp.*dJm);
e3 = c.orders.*c.beta./c.a.*Hm.*(1./c.kperpi.^2 - 1./c.kperp.^2);
e4 = -i.*c.k.*(-c.mui./c.kperpi.*dJmi./Jmi.*Hm + c.mu./c.kperp.*dHm);
e5 = -i.*c.k.*(c.epi./c.kperpi.*dJmi./Jmi.*Jm - c.ep./c.kperp.*dJm);
e6 = c.orders.*c.beta./c.a.*Jm.*(1./c.kperpi.^2 - 1./c.kperp.^2);
e7 = i.*c.k.*(c.epi./c.kperpi.*dJmi./Jmi.*Hm - c.ep./c.kperp.*dHm);
e8 = -c.orders.*c.beta./c.a.*Hm.*(1./c.kperpi.^2 - 1./c.kperp.^2);

det = e1.*e6 - e2.*e5;

Mee = (e6.*e3-e2.*e7)./det;
Meh = (e6.*e4-e2.*e8)./det;
Mhe = (e1.*e7-e3.*e5)./det;
Mhh = (e1.*e8-e4.*e5)./det;
