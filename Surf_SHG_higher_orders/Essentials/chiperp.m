function [chi2] = chiperp(lam,epifun)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
me=9.31e-31;
e=1.6e-19;
c=3e8;
omg=(2*pi/lam)*c;
chi2=-(e/(4*me*omg^2))*(epifun-1);
end

