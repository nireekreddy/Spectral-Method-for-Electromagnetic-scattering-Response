function [Er0, Et0, Ex0, Ey0] = E0cmplt_orders(x, y, c,fldfun,chi)
x(x==0 & y==0) = x(x==0 & y==0) + eps; % coordinate singularity
[th, r] = cart2pol(x, y);
a=c.a;
numfac = chi*0.25i*(pi*(a))/(c.ep);
prefac=numfac*(fldfun);

ext = r > a;
Er0=zeros(size(r));Et0=zeros(size(r));
for s = 1:length(c.orders)
expk = exp(1i.*c.orders(s).*th(ext));
expki = exp(1i.*c.orders(s).*th(~ext));
%% postive order
Er0(ext)= Er0(ext)+(prefac(s)).*(2*c.orders(s)^2./(a.*r(ext))).*(expk).*(besselj(c.orders(s),c.kperp*a).*besselh(c.orders(s),c.kperp.*r(ext)));
Et0(ext)=Et0(ext)+(prefac(s)).*(1i*c.orders(s)*c.kperp./a).*(expk).*(besselj(c.orders(s),c.kperp.*a).*(besselh(c.orders(s)-1,c.kperp.*r(ext))-besselh(c.orders(s)+1,c.kperp.*r(ext))));
% % % 
Er0(~ext)= Er0(~ext)+(prefac(s)).*(2*c.orders(s)^2./(a.*r(~ext))).*(expki).*(besselh(c.orders(s),c.kperp*a).*besselj(c.orders(s),c.kperp.*r(~ext)));
Et0(~ext)=Et0(~ext)+(prefac(s)).*(1i*c.orders(s)*c.kperp./a).*(expki).*(besselh(c.orders(s),c.kperp.*a).*(besselj(c.orders(s)-1,c.kperp.*r(~ext))-besselj(c.orders(s)+1,c.kperp.*r(~ext))));

end
Ex0 = cos(th).*Er0 - sin(th).*Et0;
Ey0 = sin(th).*Er0 + cos(th).*Et0;