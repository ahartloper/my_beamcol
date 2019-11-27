function B = b_matrix(x, L)
% Returns the B matrix evaluated at x for beam-column elements.
% @param double x: Coordinate along the length of the element
% @param double L: Length of the element
% @returns matrix: (double, 2x6) Strain-displacement matrix.
% 
% Notes:
%   - See slide 14 from week 8
B = zeros(2, 6);
nd = shape_deriv(x, L);
B(1, 1) = nd(1);
B(1, 4) = nd(2);
B(2, 2) = nd(3);
B(2, 3) = nd(4);
B(2, 5) = nd(5);
B(2, 6) = nd(6);
end
