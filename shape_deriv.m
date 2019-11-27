function Nd = shape_deriv(x, L)
% Returns the derivatives of the shape functions for beam column elements evaluated at x.
% @param double x: Coordinate along the length of the element
% @param double L: Length of the element
% @returns matrix: (double, 6x1) Dervivaties of axial, translational, rotational
%                   functions, see the Notes for the derivative orders
%
% Notes:
%   - First derivative of N1, N2 and second derivatives of N3 to N6
%   - See slide 14 from lecture 6
%   - Order of functions is: [psi_1; psi_2; phi_1; phi_2; phi_3; phi_4]
Nd = [-1. / L; 
    1. / L;
    (-6.* L + 12. *x)/L^3;
    (-4.* L + 6.* x)/L^2;
    (6.* L - 12.* x)/L^3;
    -((2* (L - 3* x))/L^2)];
end
