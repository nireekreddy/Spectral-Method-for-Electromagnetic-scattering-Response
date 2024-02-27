function value = besselh1d(nu, z)

value = 0.5.*(besselh(nu - 1, 1, z) - besselh(nu + 1, 1, z));
