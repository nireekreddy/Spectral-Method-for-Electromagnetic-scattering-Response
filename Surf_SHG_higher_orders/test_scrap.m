% [Er0, Et0] = E0cmplt_orders(x, y, sec,fldfun,chi2);
% figure;
% one=subplot(2,1,1);
% clim =(max(max(real(Er0))))*[-1 1];
% pcolor(x/1e-9, y/1e-9, real(Er0));
% shading flat;axis square;colormap(jet);colorbar
% set(gca,'fontsize',35);caxis(clim);
% two=subplot(2,1,2);
% clim =(max(max(real(Et0))))*[-1 1];
% pcolor(x/1e-9, y/1e-9, real(Et0));
% shading flat;axis square;colormap(jet);colorbar
% set(gca,'fontsize',35);caxis(clim);
% 
% %% checking the disc from Mie and SDM
% rin=a-1e-8*a; rot=a+1e-8*a; thth=pi/8;
% [xin,yin]=pol2cart(thth,rin);[xot,yot]=pol2cart(thth,rot);
% [~, Etmin, ~, ~, ~, ~] = fieldpt(xin, yin, sec);
% [~, Etmot, ~, ~, ~, ~] = fieldpt(xot, yot, sec);
% Etmot-Etmin
% [~, Et0in] = E0cmplt_orders(xin, yin, sec,fldfun,chi2);
% [~, Et0ot] = E0cmplt_orders(xot, yot, sec,fldfun,chi2);
% Et0ot-Et0in
clc;
test=sec;
test.orders=sec.orders(1);
test.AmE = []; test.AmH = [];
test.BmE = []; test.CmE = []; 
test.CmH=sec.CmH(1);test.BmH=sec.BmH(1);
epin = disproots(test, 1);
test.epi = epin;
test= modefield(test)



%%%%%
% % 
% % orders: -6
% %         ep: 1
% %         mu: 1
% %          k: 1.0472e+07
% %       beta: 0
% %        phi: []
% %          a: 1.0000e-08
% %        epi: -1.0003
% %        mui: 1
% %          x: 0
% %          y: 0
% %        AmE: 0
% %        BmE: 0
% %        CmE: 0
% %        AmH: 0
% %        BmH: -9.1986e-04 - 3.3808e-19i
% %        CmH: 2.1875e+01 - 5.9534e+16i