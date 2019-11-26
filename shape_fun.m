function N = shape_fun(x, L)
% Returns the shape functions for beam column elements evaluated at x.
% @param double x: Coordinate along the length of the element
% @param double L: Length of the element
% @returns matrix: (double, 6x1) Axial, translational, rotational functions
%
% Notes:
%   - See slide 12 from lecture 6
%   - Order of functions is: [psi_1; psi_2; phi_1; phi_2; phi_3; phi_4]
N = [1. - x/L; 
    x/L;
    1 - 3.*(x/L)^2 + 2.*(x/L)^3;
    x*(1. - x/L)^2;
    3.*(x/L)^2 - 2.*(x/L)^3;
    x*((x/L)^2 - x/L)];
end