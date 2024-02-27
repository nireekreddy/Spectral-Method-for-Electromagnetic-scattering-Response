function grand = intfunmod2(x, y, cyl1,cyl2, code)

% provides integrand using combinations of fields
% where the first field corresponds to the adjoint

% evaluate field at point
[Ex2, Ey2, Ez2] = fieldptacart(x, y, cyl2); % sec cylin adjoint
[Ex1, Ey1, Ez1] = fieldptcart(x, y, cyl1); % fun cylin

% evaluate integrand according to input
grand = eval(code);
