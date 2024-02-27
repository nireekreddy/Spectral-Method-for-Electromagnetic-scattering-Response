function [root, n] = newton(fn, cyl, init, trustr)
% parameters
tolstop = 1e-14;
tolwarn = 1e-10;
iterlim = 20; 

% initialize iteration variables
n = 0;
prev = init;

while true
  root = prev - fn(cyl, prev);

  % exit conditions
  diff = abs(root - prev)/abs(prev);
  if diff < tolstop
    break
  elseif n == iterlim;
    break
  else
    prev = root;
    n = n + 1;
  end
end

% errors and warnings
if nargin == 4 && (abs(root-init) > trustr || ~isfinite(root))
  error('Trust region violated')
elseif n == iterlim && diff > tolwarn;
  warning('Newton:IterLimitReached', ['Iteration limit was reached. Relative size of final step was ' num2str(diff) '.'])
end
