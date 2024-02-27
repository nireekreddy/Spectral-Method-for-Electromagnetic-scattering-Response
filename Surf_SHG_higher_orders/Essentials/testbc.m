function pass = testbc(orders, AmE, BmE, CmE, AmH, BmH, CmH, ep, epi, mu, mui, a, k, beta)

% Test boundary conditions using coefficients

% propagation constants
kperp = sqrt(k^2*ep*mu - beta^2);
kperpi = sqrt(k^2*epi*mui - beta^2);

nka = kperp.*a;
nika = kperpi.*a;

% plot points
theta = linspace(0, 2*pi, 10);

% output
pass = 1;

% Ez test
in = CmE.*besselj(orders, nika);
out = AmE.*besselj(orders, nka) + BmE.*besselh1(orders, nka);
pass = pass && test(in, out);

% Hz test
in = CmH.*besselj(orders, nika);
out = AmH.*besselj(orders, nka) + BmH.*besselh1(orders, nka);
pass = pass && test(in, out);

% Etheta test
in = i./kperpi.^2.*(beta./a.*i.*orders.*CmE.*besselj(orders, nika) - k.*mui.*kperpi.*CmH.*besseljd(orders, nika));
out = i./kperp.^2.*(beta./a.*i.*orders.*(AmE.*besselj(orders, nka) + BmE.*besselh1(orders, nka)) - k.*mu.*kperp.*(AmH.*besseljd(orders, nka) + BmH.*besselh1d(orders, nka)));
pass = pass && test(in, out);

% Htheta test
in = i./kperpi.^2.*(beta./a.*i.*orders.*CmH.*besselj(orders, nika) + k.*epi.*kperpi.*CmE.*besseljd(orders, nika));
out = i./kperp.^2.*(beta./a.*i.*orders.*(AmH.*besselj(orders, nka) + BmH.*besselh1(orders, nka)) + k.*ep.*kperp.*(AmE.*besseljd(orders, nka) + BmE.*besselh1d(orders, nka)));
pass = pass && test(in, out);

% Er test
in = i./kperpi.^2.*(beta.*kperpi.*CmE.*besseljd(orders, nika) + k.*mui./a.*i.*orders.*CmH.*besselj(orders, nika));
out = i./kperp.^2.*(beta.*kperp.*(AmE.*besseljd(orders, nka) + BmE.*besselh1d(orders, nka)) + k.*mu./a.*i.*orders.*(AmH.*besselj(orders, nka) + BmH.*besselh(orders, nka)));
pass = pass && test(epi.*in, ep.*out);

% Hr test
in = i./kperpi.^2.*(beta.*kperpi.*CmH.*besseljd(orders, nika) - k.*epi./a.*i.*orders.*CmE.*besselj(orders, nika));
out = i./kperp.^2.*(beta.*kperp.*(AmH.*besseljd(orders, nka) + BmH.*besselh1d(orders, nka)) - k.*ep./a.*i.*orders.*(AmE.*besselj(orders, nka) + BmE.*besselh(orders, nka)));
pass = pass && test(mui.*in, mu.*out);

function pass = test(x, y)
tol = eps(single(1));
res = norm(sum(x-y));
pass = res < tol;
