function grand = intfunmod(x, y, cyl,cyl0,fldfun,chi,code)

% provides integrand using combinations of fields
% where the first field corresponds to the adjoint

% evaluate field at point
[Ex1, Ey1, ~] = fieldptacart(x, y, cyl);
[~, ~,Ex2, Ey2] = E0cmplt_orders(x, y, cyl0,fldfun,chi); % E0

% evaluate integrand according to input
grand = eval(code);
