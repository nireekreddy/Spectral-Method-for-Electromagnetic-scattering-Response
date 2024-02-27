function [vals, err] = contour(cyl, c, r, o, phi, nroots, zsings, zroots, tol)

% initial number of points
pts = 32;
% maximum number of points
maxpts = 2^14;

% pre-allocate variables
[t, expt, dexpt, z, fnval] = deal(zeros([1 2^11]));

% set up evaluation points
t(1:pts+1) = linspace(0,1,pts+1);
t(pts+1) = 0;

% define additional points
add = 1:pts;

% count number of singularities and known roots in contour
nsings = length(zsings); kroots = length(zroots);

% loop until count of roots reaches tolerance
while true
  % create oblate contour
  expt(add) = exp(1i*phi).*(cos(2.*pi.*t(add)) + 1i.*o.*sin(2.*pi.*t(add)));
  dexpt(add) = exp(1i.*phi).*(o.*cos(2.*pi.*t(add)) + 1i.*sin(2.*pi.*t(add)));
  z(add) = c + r.*expt(add);
  
  % evaluate ratio of function and derivative
  fnval(add) = dexpt(add)./cyldispnewt(cyl, z(add));
  
  % evaluate error based on count of roots
  cnt = r/pts*sum(fnval(1:pts));
  err = abs(cnt - nroots + nsings - kroots);

  % end loop conditions
  if err < tol
    break
  elseif pts > maxpts
    warning(['Maximum number of intgration points, ' num2str(maxpts) ', exceeded. Root count is ' num2str(sign(real(cnt))*abs(cnt)+nsings-kroots) '.'])
    break
  else % intercalate points
    pts = 2*pts;
    add = pts/2+1:pts;
    t(add) = t(1:pts/2)+1/pts;
  end
end

% evaluate the sum of powers
s = 1:nroots;
% sp = r/pts*sum(expt(1:pts).^s.*fnval(1:pts));
sp = r./pts.*sum(bsxfun(@times, bsxfun(@power, expt(1:pts).', s), fnval(1:pts).'));

% sum of powers of singularities
if nsings > 0
  % normalized coordinates
  us = (zsings-c)./r;
  % sps = 2*sum(us.^s);
  sps = sum(bsxfun(@power, us.', s), 1);
  
  % deflate singularities
  sp = sp + sps;
end

% sum of powers of known roots
if kroots > 0
  % normalized coordinates
  ur = (zroots-c)./r;
  % spr = 2*sum(us.^s);
  spr = sum(bsxfun(@power, ur.', s), 1);
  
  % deflate roots
  sp = sp - spr;
end

% recursive evaluation symmetric polynomials
% index begins with zeroth order
ep = zeros([1 nroots]);
ep(1) = 1;
for k = 1:nroots
  inds = 1:k;
  ep(k+1) = 1/k * sum((-1).^(inds-1).*ep(k-inds+1).*sp(inds));
end

% define polynomial and find roots
poly = (-1).^(0:nroots) .* ep;
uroots = roots(poly);

% shift roots from unit circle
vals = c + r.*uroots;
