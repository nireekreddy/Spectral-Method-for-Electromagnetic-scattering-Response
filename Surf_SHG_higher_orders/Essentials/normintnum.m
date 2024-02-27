function fac = normintnum(cyl) 

%function fac = normintnum(cyl,cyl2,ffstr)

% define circular integration domain
xmax = cyl.a; xmin = -cyl.a;
ymax = @(x) sqrt(cyl.a^2-x.^2);
ymin = @(x) -sqrt(cyl.a^2-x.^2);

% intgrand
EtEt = @(x, y) intfun(x, y, cyl, cyl, 'Ex1.*Ex2+Ey1.*Ey2+Ez1.*Ez2');
%EtEt = @(x, y) intfunc(x, y, cyl, cyl, 'conj(Ex1).*Ex2+conj(Ey1).*Ey2+conj(Ez1).*Ez2');

% integration factor
intotal = integral2(EtEt,xmin,xmax,ymin,ymax);
fac = sqrt(intotal);
