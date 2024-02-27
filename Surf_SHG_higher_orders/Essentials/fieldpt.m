function [Er, Et, Ez, Hr, Ht, Hz] = fieldpt(x, y, c)

% Obtains Cartesian fields components at specified point

% plot points
x(x==0 & y==0) = x(x==0 & y==0) + eps; % coordinate singularity
[th, r] = cart2pol(x-c.x, y-c.y);

% test if plot points exterior to cylinder
ext = r > c.a;
c = c.initcoeff;

[Ez, Er, Et, Hz, Hr, Ht] = deal(zeros(size(r)));
for s = 1:length(c.orders)
  bessj = besselj(c.orders(s),c.kperp.*r(ext));
  bessh = besselh1(c.orders(s),c.kperp.*r(ext));
  bessjd = besseljd(c.orders(s),c.kperp.*r(ext));
  besshd = besselh1d(c.orders(s),c.kperp.*r(ext));
  bessji = besselj(c.orders(s),c.kperpi.*r(~ext));
  bessjid = besseljd(c.orders(s),c.kperpi.*r(~ext));

  expk = exp(i.*c.orders(s).*th(ext));
  expki = exp(i.*c.orders(s).*th(~ext));

  Ez(ext) = Ez(ext) + (c.AmE(s).*bessj+c.BmE(s).*bessh).*expk;
  Et(ext) = Et(ext) + i./c.kperp.^2.*(c.beta./r(ext).*i.*c.orders(s).*(c.AmE(s).*bessj+c.BmE(s).*bessh)-c.k.*c.mu.*c.kperp.*(c.AmH(s).*bessjd+c.BmH(s).*besshd)).*expk;
  Er(ext) = Er(ext) + i./c.kperp.^2.*(c.beta.*c.kperp.*(c.AmE(s).*bessjd+c.BmE(s).*besshd)+c.k.*c.mu./r(ext).*i.*c.orders(s).*(c.AmH(s).*bessj+c.BmH(s).*bessh)).*expk;

  Ez(~ext) = Ez(~ext) + c.CmE(s).*bessji.*expki;
  Et(~ext) = Et(~ext) + i./c.kperpi.^2.*(c.beta./r(~ext).*i.*c.orders(s).*c.CmE(s).*bessji - c.k.*c.mui.*c.kperpi.*c.CmH(s).*bessjid).*expki;
  Er(~ext) = Er(~ext) + i./c.kperpi.^2.*(c.beta.*c.kperpi.*c.CmE(s).*bessjid+c.k.*c.mui./r(~ext).*i.*c.orders(s).*c.CmH(s).*bessji).*expki;

  Hz(ext) = Hz(ext) + (c.AmH(s).*bessj+c.BmH(s).*bessh).*expk;
  Ht(ext) = Ht(ext) + i./c.kperp.^2.*(c.beta./r(ext).*i.*c.orders(s).*(c.AmH(s).*bessj+c.BmH(s).*bessh)+c.k.*c.ep.*c.kperp.*(c.AmE(s).*bessjd+c.BmE(s).*besshd)).*expk;
  Hr(ext) = Hr(ext) + i./c.kperp.^2.*(c.beta.*c.kperp.*(c.AmH(s).*bessjd+c.BmH(s).*besshd)-c.k.*c.ep./r(ext).*i.*c.orders(s).*(c.AmE(s).*bessj+c.BmE(s).*bessh)).*expk;

  Hz(~ext) = Hz(~ext) + c.CmH(s).*bessji.*expki;
  Ht(~ext) = Ht(~ext) + i./c.kperpi.^2.*(c.beta./r(~ext).*i.*c.orders(s).*c.CmH(s).*bessji + c.k.*c.epi.*c.kperpi.*c.CmE(s).*bessjid).*expki;
  Hr(~ext) = Hr(~ext) + i./c.kperpi.^2.*(c.beta.*c.kperpi.*c.CmH(s).*bessjid-c.k.*c.epi./r(~ext).*i.*c.orders(s).*c.CmE(s).*bessji).*expki;
end

% Cartesian fields
%Ex = cos(th).*Er - sin(th).*Et;
%Ey = sin(th).*Er + cos(th).*Et;
