function fs = assemble_fs(sigma, A, YZ)
% Returns the section force vector at a given x.
% @param matrix sigma: (double, nx1) Stress for each fiber at x
% @param matrix A: (double, nx1) Area of each fiber at x
% @param matrix YZ: (double, nx2) y,z coordinates of each fiber
% @returns matrix: (double, 2x1) Resisting force vector for the cross-section at x
%
% Notes:
%   - n is the number of fibers for the cross-section.
%   - This function should be evaluated at a know value of x, the input
%   parameters should correspond to this value of x.
%   - ks corresponds to: [P, M] = ks * [\epsilon, \phi]
n = length(sigma);
fs = zeros(2, 1);
for i = 1:n
    l = l_vec(YZ(i, :));
    fs = fs + sigma(i) * A(i) * l';
end
end
