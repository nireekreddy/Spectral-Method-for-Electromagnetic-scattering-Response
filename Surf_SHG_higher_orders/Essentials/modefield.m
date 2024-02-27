function cyl = modefield(cyl)
% determines and normalizes fields of a mode

% set zero incidence
cyl.AmE = 0;
cyl.AmH = 0;

% test for TM or TE mode
if cyl.orders == 0 || cyl.beta == 0
  [~, TEzero, TMzero] = cyldisp(cyl);
end  

% Bessel arguments
ka = cyl.kperp*cyl.a;
kia = cyl.kperpi*cyl.a;

% get modal fields
if cyl.orders ~= 0 && cyl.beta ~= 0
  % hybrid mode for m ~= 0
  cyl.CmE = 1;
  cyl.CmH = i*cyl.k*cyl.a/cyl.orders/cyl.beta * (cyl.ep/cyl.kperp*besselh1d(cyl.orders, ka)/besselh(cyl.orders, ka) - cyl.epi/cyl.kperpi*besseljd(cyl.orders, kia)/besselj(cyl.orders, kia)) / (cyl.kperp^-2 - cyl.kperpi^-2);
elseif abs(TMzero) < abs(TEzero)
  % enforces Ez mode for m = 0
  cyl.CmE = 1; cyl.CmH = 0;
else
  % enforces Hz mode for m = 0
  cyl.CmE = 0; cyl.CmH = 1;
end

% normalization
fac = normint(cyl);
cyl.CmE = cyl.CmE/fac; cyl.CmH = cyl.CmH/fac;

% exterior fields
cyl.BmE = besselj(cyl.orders, cyl.kperpi*cyl.a)/besselh(cyl.orders, cyl.kperp*cyl.a) * cyl.CmE;
cyl.BmH = besselj(cyl.orders, cyl.kperpi*cyl.a)/besselh(cyl.orders, cyl.kperp*cyl.a) * cyl.CmH;
end
