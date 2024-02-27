clear; clf;clc
cyl = CylinderGeometry;

cyl.a =30*1e-9;
lam=400*1e-9;


% simulation properties
cyl.orders = -4:4;

% background properties
cyl.ep = 1; cyl.mu = 1.0;

% field properties
cyl.k = 2*pi/lam; cyl.beta = 0.0;
cyl.phi = 0;

% cylinder properties
cyl.epi = drude(lam);
cyl.mui = 1;

% incident field Bessel coefficients
E0 = 0; H0 = 1;
cyl.AmE = E0.*i.^cyl(1).orders.*exp(-i.*cyl(1).orders.*cyl.phi); 
cyl.AmH = H0.*i.^cyl(1).orders.*exp(-i.*cyl(1).orders.*cyl.phi); 

% inverse Mie coefficients
[Mee, Meh, Mhe, Mhh] = miecoeff(cyl);

% Mie field coefficients
B = [diag(Mee) diag(Meh); diag(Mhe) diag(Mhh)]\[cyl.AmE.'; cyl.AmH.'];
cyl.BmE = B(1:length(B)/2).'; cyl.BmH = B(length(B)/2+1:end).';
cyl.CmE = intcoeff(cyl.orders, cyl.AmE, cyl.BmE, cyl.kperp, cyl.kperpi, cyl.a);
cyl.CmH = intcoeff(cyl.orders, cyl.AmH, cyl.BmH, cyl.kperp, cyl.kperpi, cyl.a);

cyl.AmH=[];
% calculate total fields
[x, y] = meshgrid(linspace(-5*cyl.a,5*cyl.a,800), linspace(-5*cyl.a,5*cyl.a,800));
[Er, Et, Ez, Hr, Ht, Hz] = fieldpt(x, y, cyl);

%%%%cartesian
[th, r] = cart2pol(x, y);
Ex = cos(th).*Er - sin(th).*Et;
Ey = sin(th).*Er + cos(th).*Et;
Hx = cos(th).*Hr - sin(th).*Ht;
Hy = sin(th).*Hr + cos(th).*Ht;
Z=377;
%%%plot
figure;
clim = max(max(real(Ez)))*[-1 1];
subplot(3,2,1)
pcolor(x, y, real(Ez))
axis square
shading flat
colorbar
caxis(clim);
subplot(3,2,2)
pcolor(x, y, real(Hz)/Z)
axis square
shading flat
colorbar
caxis(clim);
subplot(3,2,3)
pcolor(x, y, real(Er))
axis square
shading flat
colorbar
caxis(clim);
subplot(3,2,4)
pcolor(x, y, real(Hx)/Z)
axis square
shading flat
colorbar
caxis(clim);
subplot(3,2,5)
pcolor(x, y, real(Et))
axis square
shading flat
colorbar
caxis(clim);
subplot(3,2,6)
pcolor(x, y, real(Hy)/Z)
axis square
shading flat
colorbar
caxis(clim);


%% plot vrosection
% % figure;
% % yy=linspace(2*cyl.a,-2*cyl.a,1001);
% % xx=zeros(1,length(yy));
% % [Erc, Etc, Ezc, Hrc, Htc, Hzc] = fieldpt(xx, yy, cyl);
% % Z=377;
% % plot(yy/1e-9,real(Erc));