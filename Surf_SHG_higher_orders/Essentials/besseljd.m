function value = besseljd(nu, z)

value = 0.5.*(besselj(nu - 1, z) - besselj(nu + 1, z));
