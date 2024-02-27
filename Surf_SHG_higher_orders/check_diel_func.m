clc;clear;lamvec=[];
lami=420e-9;lamf=1600e-9;nfac=600;dlam=(lamf-lami)/nfac;epigold=[];epidr=[];
for uu=1:nfac;
lam=lami+(uu-1)*dlam % free space wavelength
lamvec=[lamvec lam];
lammic=lam/1e-6;
kfun = 2*pi/lam; ksec = 2*kfun; % free space wavelengths at FF and SH    
epigold = [epigold au(lammic)]; 
epidr=[epidr drude(lam)];
end
subplot(211)
plot(lamvec/1e-9,real(epigold),'r',lamvec/1e-9,real(epidr),'b');
subplot(212)
plot(lamvec/1e-9,imag(epigold),'r',lamvec/1e-9,imag(epidr),'b')