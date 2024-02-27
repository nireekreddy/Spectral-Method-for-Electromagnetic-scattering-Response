figure;
one=subplot(2,2,1);
pcolor(x/1e-9, y/1e-9, real(Drs)./(max(max(real(Drs)))));
set(gca,'TickLabelInterpreter','latex')
shading flat;axis square;colormap(jet);colorbar
set(gca,'fontsize',35);
hold on
two=subplot(2,2,3);
pcolor(x/1e-9, y/1e-9, real(Ets)/(max(max(real(Ets)))));
grid on;shading flat;axis square;colorbar
colormap(jet);set(gca,'fontsize',35)
set(gca,'TickLabelInterpreter','latex')

three=subplot(2,2,2);
%clim =(max(max(abs(Drsberg-Drs))))*[0 1];
pcolor(x/1e-9, y/1e-9,abs(Drsberg-Drs)./abs(Drs));
shading flat;axis square;
set(gca,'TickLabelInterpreter','latex')
colorbar
colormap(jet)
%caxis(clim);
set(gca,'fontsize',35)
box('on')
%%%%%%%%%%%%%%%
four=subplot(2,2,4)
%clim =(max(max(abs(Etsberg-Ets))))*[0 1];
pcolor(x/1e-9, y/1e-9, abs(Etsberg-Ets)./abs(Ets));
set(gca,'TickLabelInterpreter','latex')
shading flat;axis square;
colorbar
colormap(jet)
%caxis(clim);
set(gca,'fontsize',35)
box('on')
