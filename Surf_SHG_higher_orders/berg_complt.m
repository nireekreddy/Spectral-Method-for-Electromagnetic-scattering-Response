% number of radial eigenmodes
%% evaluating E0 for botht the orders
[Er0, Et0,~,~] = E0cmplt_orders(x, y, sec,fldfun,chi2);
tmode = sec;
%% modes and overlap
[Ern, Etn] = deal(zeros(size(x)));
for ll=1:length(sec.orders)
 tmode.orders=sec.orders(ll);
tmode.AmE = []; tmode.AmH = [];
tmode.BmE = []; tmode.CmE = []; 
tmode.CmH=sec.CmH(ll);tmode.BmH=sec.BmH(ll);
% dispersion relation root search and modified cylinder
epin = disproots(tmode, 1);
tmode.epi = epin;
tmode= modefield(tmode);
% overlap for each radial mode
olap=normintnummod(tmode,sec,fldfun,chi2); %integrating numerically
detn=(episec-tmode.ep)./((tmode.epi-episec));
faccrnt=olap.*detn;
% % Adding all the radial modes
[Er, Et, ~, ~, ~, ~] = fieldpt(x, y, tmode);
                 Ern = Ern +faccrnt.*Er; 
                 Etn = Etn +faccrnt.*Et; 
end

% total second harmonic field - from modes and E0
Ersberg =(Er0+Ern); Etsberg =(Et0+Etn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%calculating the radial displacement field
Drsberg = Ersberg;
Drsberg(r>sec.a) = Drsberg(r>sec.a).*ep;
Drsberg(r<=sec.a) = Drsberg(r<=sec.a).*episec;

%%%%%%Plotting the fields
hold on
subplot(2,3,2)
clim =(max(max(real(Drsberg))))*[-1 1];
pcolor(x, y,real(Drsberg));
title('SDM-$D_{\rho}$','Interpreter','latex')
shading flat;axis square;
colorbar
colormap(jet)
caxis(clim);
set(gca,'XTick',[],'YTick',[]);
set(gca,'fontsize',25)
box('on')
%%%%%%%%%%%%%%%
subplot(2,3,5)
clim =(max(max(real(Etsberg))))*[-1 1];
pcolor(x, y, real(Etsberg));
title('SDM-$E_{\theta}$','Interpreter','latex')
shading flat;axis square;
colorbar
colormap(jet)
caxis(clim);
set(gca,'XTick',[],'YTick',[]);
set(gca,'fontsize',25)
%% plotting the difference
hold on
subplot(2,3,3)
clim =(max(max(abs(Drsberg-Drs))))*[0 1];
pcolor(x, y,abs(Drsberg-Drs));
title('SDM-$D_{\rho}$','Interpreter','latex')
shading flat;axis square;
colorbar
colormap(jet)
caxis(clim);
set(gca,'XTick',[],'YTick',[]);
set(gca,'fontsize',25)
box('on')
%%%%%%%%%%%%%%%
subplot(2,3,6)
clim =(max(max(abs(Etsberg-Ets))))*[0 1];
pcolor(x, y, abs(Etsberg-Ets));
title('SDM-$E_{\theta}$','Interpreter','latex')
shading flat;axis square;
colorbar
colormap(jet)
caxis(clim);
set(gca,'XTick',[],'YTick',[]);
set(gca,'fontsize',25)

