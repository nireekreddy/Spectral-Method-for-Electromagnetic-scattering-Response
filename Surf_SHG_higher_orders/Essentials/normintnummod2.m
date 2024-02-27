function intotal = normintnummod2(cyl1,cyl2) 

%  cyl1 fun cylinder 

% define circular integration domain
xmax = cyl1.a; xmin = -cyl1.a;
ymax = @(x) sqrt(cyl1.a^2-x.^2);
ymin = @(x) -sqrt(cyl1.a^2-x.^2);

% intgrand
EtEt = @(x, y) intfunmod2(x, y, cyl1,cyl2,'(Ex1).^2.*Ex2+(Ey1).^2.*Ey2+(Ez1).^2.*Ez2');
%EtEt = @(x, y) intfunc(x, y, cyl, cyl, 'conj(Ex1).*Ex2+conj(Ey1).*Ey2+conj(Ez1).*Ez2');
% integration factor
intotal = integral2(EtEt,xmin,xmax,ymin,ymax);
