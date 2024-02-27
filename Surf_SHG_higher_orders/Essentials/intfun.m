function grand = intfun(x, y, cyl1, cyl2, code)

% provides integrand using combinations of fields
% where the first field corresponds to the adjoint

% evaluate field at point
[Ex1, Ey1, Ez1] = fieldptacart(x, y, cyl1);
[Ex2, Ey2, Ez2] = fieldptcart(x, y, cyl2); % E0

% evaluate integrand according to input
grand = eval(code);
