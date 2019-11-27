function ad = ad_matrix(x, L)
% Returns the a_d matrix evaluated at x for beam-column elements.
% @param double x: Coordinate along the length of the element
% @param double L: Length of the element
% @returns matrix: (double, 2x6) Cross-section displacement matrix.
% 
% Notes:
%   - See slide 11 from week 8
ad = zeros(2, 6);
n = shape_fun(x, L);
ad(1, 1) = n(1);
ad(1, 4) = n(2);
ad(2, 2) = n(3);
ad(2, 3) = n(4);
ad(2, 5) = n(5);
ad(2, 6) = n(6);
end
