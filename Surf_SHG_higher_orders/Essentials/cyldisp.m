function [detm, fac1, fac2] = cyldisp(cyl, epi)

% set epi if not provided
if nargin == 1
  epi = cyl.epi;
end

% interior exterior arguments
kperp = sqrt(cyl.k.^2.*cyl.ep.*cyl.mu - cyl.beta.^2);
kperpi = sqrt(cyl.k.^2.*epi.*cyl.mu - cyl.beta.^2);
ka = kperp.*cyl.a;
kia = kperpi.*cyl.a;

% equation ratios
jrat = besseljd(cyl.orders, kia)./besselj(cyl.orders, kia)./kia;
hrat = besselh1d(cyl.orders, ka)./besselh(cyl.orders, ka)./ka;

detm = cyl.k.^2.*(jrat-hrat).*(epi.*jrat-cyl.ep.*hrat) - cyl.orders.^2*cyl.beta.^2.*(kia.^-2-ka.^-2).^2;

% individual factors to determine TM/TE mode
fac1 = jrat-hrat;
fac2 = epi.*jrat-cyl.ep.*hrat;
