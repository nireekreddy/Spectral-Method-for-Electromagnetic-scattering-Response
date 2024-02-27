function [ epsmet ] = drude(lam)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
omg= 1239.83/(lam/1e-9);
omgp = 9.1 ; %eV
gam=0.072; %eV damping
epsinf=9;
epsmet=epsinf-omgp^2/(omg^2+i*omg*gam);
end