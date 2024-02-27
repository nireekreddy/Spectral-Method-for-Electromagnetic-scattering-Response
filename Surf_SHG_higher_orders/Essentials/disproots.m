function roots = disproots(cyl, nroots)

% number of contours
% each singularity known to have 2 roots nearby
sings = ceil(nroots/2);

% get zero of J, corresponds to kperpi*a
jzero = besselzero(cyl.orders, max(4, sings+1));
% insert ficticious singularity for delineation purposes
ising = interp1(1:4, jzero(1:4), 0, 'pchip');
% calculates corresponding epsilon
epis = (([ising jzero]./cyl.a).^2 + cyl.beta.^2)./cyl.k.^2./cyl.mu;

% find plasmonic root if it exists
if cyl.orders == 0 && real(cyl.kperp) ~= 0
  proot = [];
else
  proot = newton(cyl, -1-1e-7i, 1e6);
end

% boundary is 1/4 distance between singularities
dz = diff(epis);
llim = epis(1:end-1)+dz./4;

% extend boundary if lowest boundary if close to unity
if llim(1) < 50; 
  if length(proot) == 0 || abs(proot) > 20;
    % contour excludes plasmonic root
    llim(1) = -10;
  else
    llim(1) = -abs(proot)-10;
  end
end

% calculate center and radius of contours
rad = diff(llim)./2;
cen = llim(1:end-1) + rad;

droots = zeros([1 2*sings]);
ncntrs = ceil((nroots-length(proot))/2);
% initiate root statistics based on first two pairs of roots
for k = 1:min(2,ncntrs)
  droots([2*k-1 2*k]) = search(cyl, epis(k+1), proot, cen(k), rad(k));
end

% targeted search based on root statistics
for k = 3:ncntrs
  droots([2*k-1 2*k]) = predict(cyl, epis(k-1:k+1), droots(2*k-5:2*k-2));
end

% combine roots and exclude unwanted roots
roots = [proot droots];
roots = roots(1:nroots);

function root = newton(cyl, init, trustr)
% parameters
tol = 1e-14;
stol = 1e-10;
iterlim = 20; 

% initialize iteration variables
n = 0;
prev = init;

while true
  root = prev - cyldispnewt(cyl, prev);

  % exit conditions
  diff = abs(root - prev)/abs(prev);
  if diff < tol
    break
  elseif n > iterlim;
    break
  else
    prev = root;
    n = n + 1;
  end
end

% errors and warnings
if abs(root-init) > trustr || ~isfinite(root)
  error('Trust region violated')
elseif n > iterlim
  warning('Iteration limit was reached.')
end

function roots = search(cyl, sing, proot, cen, rad)
% specify parameters
crstol = 1e-5;
fintol = 1e-10;
consize = 1.5;
maxobl = 1/10;
trustr = 10;

% include polynomial singularity if it exists in large contour
% bessel singularities known to be double poles
if (cyl.orders ~= 0 && cyl.beta ~= 0) && abs(cen - cyl.kperp) < rad
  sings = [cyl.beta^2/cyl.k^2 sing sing];
else
  sings = [sing sing];
end

% exclude plasmonic root if not in large contour
if length(proot) ~= 0 && abs(cen - proot) > rad
  proot = [];
end

% large contour to include both roots
phi = 0; obl = 1;
roots = contour(cyl, cen, rad, obl, phi, 2, sings, proot, crstol);

% sort roots: second assumed to be near singularity
roots = sort(roots);

% polish roots using Newton-Raphson
roots(1) = newton(cyl, roots(1), trustr);

try
  % inaccuracy relative to singularity causes failure to converge
  roots(2) = newton(cyl, roots(2), trustr);
catch
  % fine contour centered near singularity to refine initial guess
  fcen = (roots(2)+sing)/2;
  frad = max(consize*abs(sing-fcen), 1e-3*rad);
  fphi = atan2(imag(sing-fcen), real(sing-fcen));
  fobl = max(consize*abs(imag(sing-fcen))/frad, maxobl);
  fsings = [sing sing]; fzeros = [];
  roots(2) = contour(cyl, fcen, frad, fobl, fphi, 1, fsings, fzeros, fintol);
  roots(2) = newton(cyl, roots(2), trustr);
end

function roots = predict(cyl, sings, proots)
trustr = 10;

roots = zeros([2 1]);

% predict location of flat root
offset = proots([1 3]) - sings(1:2);
poffset = 2*offset(2)-offset(1);
roots(1) = newton(cyl, sings(3) + poffset, trustr);

% predict location of singularity root
offset = proots([2 4]) - sings(1:2);
poffset = 2*offset(2)-offset(1);
roots(2) = newton(cyl, sings(3) + poffset, trustr);

%fdifference = roots(1) - sings(3) - poffset
%sdifference = roots(2) - sings(3) - poffset

% defunct function
function roots = direct(cyl, sing)
% specify parameters
fintol = 1e-10;

% known to be double poles
sings = kron(sing, [1 1]);

% find flat root 
%roots(1) = newton(cyl, sing/2+1e-6i);
roots(1) = newton(cyl, 1+1e-6i);

% fine contour centered near singularity
fcen = sing; frad = abs(sing-roots(1))/2;
fphi = 0; fobl = 1;
roots(2) = contour(cyl, fcen, frad, fobl, fphi, 1, sings, fintol);

% polish root using Newton-Raphson
roots(2) = newton(cyl, roots(2));

% defunct function
function roots = focussearch(cyl, sings, proots)
% specify parameters
crstol = 1e-5;
fintol = 1e-10;
consize = 1.5;
maxobl = 1/10;

% predict location of flat root
offset = proots([1 3]) - sings(1:2);
poffset = 2*offset(2)-offset(1);

% contour parameters
cen = sings(3) + poffset;
rad = 0.1*(sings(3) - sings(2));
obl = max(consize*abs(imag(sings(3)-cen))/rad, maxobl);
phi = 0;
roots(1) = contour(cyl, cen, rad, obl, phi, 1, [], crstol);

% predict location of singularity root
offset = proots([2 4]) - sings(1:2);
poffset = 2*offset(2)-offset(1);

% contour parameters
cen = sings(3) + poffset/2;
sing = kron(sings(3), [1 1]);
rad = consize*abs(poffset/2);
obl = max(consize*abs(imag(sings(3)-cen))/rad, maxobl);
phi = atan2(imag(poffset), real(poffset));
roots(2) = contour(cyl, cen, rad, obl, phi, 1, sing, fintol);

%fdifference = roots(1) - cen
%sdifference = roots(2) - sings(3) - poffset
