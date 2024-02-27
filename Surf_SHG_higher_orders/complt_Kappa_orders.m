%% Surface SHG Mie code
%%% please see last part of the code plot SH fields along the 
clear;clc;
addpath('C:\Users\kothakap\Desktop\Surf_SHG_higher_orders\Essentials')
% background properties
ep = 1; mu = 1;
%% incident field properties --  [0,exp(-i*kfun*x),0] only y component

phi = pi/2; % angle to be used in Jacobi angle expansion
beta = 0.0; % normal incidence \beta=0
%% cylinder properties at FF and SH
a = 30e-9; % cylinder radius
%fun.a=a;sec.a=a;
mui = 1; % non-magnetic

chi2=1e-20;  % switch to "chiperp(lam,epifun);" %% chi^(2)_{\prep,\prep,\perp} expression from hydrodynamic model 
%% FF cyl incident and scattered fields
orderfun = -4:4; %angular orders to be included in FF scattering
ordersec = -4:4; %% angular orders to be included in SH calc
% define incident parameters for fundamental
Hz = 1; % TE pol
lamvec=[];Delta=[];
lami=950e-9;lamf=1350e-9;nfac=1000;dlam=(lamf-lami)/nfac;
for uu=1:nfac;
lam=lami+(uu-1)*dlam % free space wavelength
lamvec=[lamvec lam];
lammic=lam/1e-6;
kfun = 2*pi/lam; ksec = 2*kfun; % free space wavelengths at FF and SH    
epifun =au(lammic); %\epsilon(\omega) of cylinder
episec =au(lammic/2); %\epsilon(2\omega) of cylinder
% define cylinder representing fundamental field
fun = CylinderGeometry;
fun.orders = orderfun;fun.k = kfun; fun.beta = beta; fun.phi = phi;
fun.ep = ep; fun.epi = epifun;fun.mu = mu; fun.mui = mui;fun.a = a;

% Jacobi-Anger
fun.AmH = Hz.*i.^fun.orders.*exp(-i.*fun.orders.*phi); 

% inverse Mie coefficients
[Mee, Meh, Mhe, Mhh] = miecoeff(fun);
% obtain scattered and interior field
fun.BmH = fun.AmH./Mhh;
fun.CmH = intcoeff(fun.orders, fun.AmH, fun.BmH, fun.kperp, fun.kperpi, fun.a);
% zero TM components
fun.AmE = []; fun.CmE = []; fun.BmE = [];

%% SH calc begins here
% define cylinder representing second harmonic field
sec = CylinderGeometry;
sec.orders = ordersec;
sec.k = ksec; sec.beta = beta;sec.ep = ep; sec.epi = episec;
sec.mu = mu; sec.mui = mui;sec.a = a;

% TM components and incident fields for SH are zero 
sec.AmE = []; sec.AmH = [];sec.BmE = []; sec.CmE = []; 
%% boundary conditions - Heinz for SH
% defining dummy cyl for FF to separate anglr ordrs
funt=CylinderGeometry; fn=[];
funt.orders = 0;funt.k = kfun; funt.beta = beta; funt.phi = phi;
funt.ep = ep; funt.epi = epifun;funt.mu = mu; funt.mui = mui;funt.a = a;

for rr = 1:length(orderfun)
% def tmp FF cyl to seperate orders
funt.orders=fun.orders(rr);funt.AmE=0;funt.BmE=0;funt.CmE=0;
funt.AmH=fun.AmH(rr);funt.BmH=fun.BmH(rr);funt.CmH=fun.CmH(rr);
% extracting FF anglr ordr fields       
[ffstr, ~, ~, ~, ~, ~]=fieldpt(a-1e-8*a, 0, funt);
fn=[fn ffstr];
end
rw=(fn);clm=wrev(fn).';
matr=clm*rw;
ms=min(sec.orders);
%% funfld for the various higher SH anglr modes
for rr = 1:length(orderfun)
fldfun(rr)=sum(diag(matr,ms+rr-1)); 
end

detn=[];
for ll = 1:length(ordersec)
% defining Hnz cndts for SH CmH and BmH
trm1 = -(sec.orders(ll)*i/(a*sec.ep))*(chi2)*(fldfun(ll));
trm2 = besselh1d(sec.orders(ll), sec.kperp*sec.a)/(sec.kperp);
trm3 = ((besselh(sec.orders(ll), sec.kperp*sec.a))/(besselj(sec.orders(ll), sec.kperpi*sec.a)))*(besseljd(sec.orders(ll), sec.kperpi*sec.a)/(sec.kperpi));
trm4 = (trm2-trm3)*(-1i*sec.k*sec.mu);
sec.BmH(ll)=trm4^(-1)*trm1; 
sec.CmH(ll)=sec.BmH(ll)*(besselh(sec.orders(ll), sec.kperp*sec.a)/besselj(sec.orders(ll), sec.kperpi*sec.a));
% % negmode = sec; negmode.orders=sec.orders(ll);
% % negmode.AmE = []; negmode.AmH = [];
% % negmode.BmE = []; negmode.CmE = []; 
% % negmode.CmH=sec.CmH(ll);negmode.BmH=sec.BmH(ll);
% % % dispersion relation root search and modified cylinder
% % epin = disproots(negmode,1);
% % detn =[detn abs(epin-episec)];
end
%% higher order angluar modes
thth=linspace(0+0.0001,2*pi-0.0001,500);rm=a-0.000000001*a;
[xm,ym]=pol2cart(thth, rm);rp=a+0.000000001*a;
[xp,yp]=pol2cart(thth, rp);
[~, Etsdm, ~, ~, ~, ~] = fieldpt(xm,ym, sec);
[~, Etsdp, ~, ~, ~, ~] = fieldpt(xp,yp, sec);
Etd=Etsdp-Etsdm;Etp=max(Etsdp+Etsdm);
Delta= [Delta max(abs(Etd./Etp))];

end
subplot(313)
plot(lamvec/1e-9,Delta*100,'k')
set(gca,'fontsize',32)
set(gca,'TickLabelInterpreter','latex')
min(Delta*100)
