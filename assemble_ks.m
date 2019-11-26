function ks = assemble_ks(E, A, YZ)
% Returns the section stiffness matrix at a given x.
% @param matrix E: (double, nx1) Tangent modulus for each fiber at x
% @param matrix A: (double, nx1) Area of each fiber at x
% @param matrix YZ: (double, nx2) y,z coordinates of each fiber
% @returns matrix: (double, 2x2) Tagent section stiffness matrix for the cross-section at x
%
% Notes:
%   - n is the number of fibers for the cross-section.
%   - This function should be evaluated at a know value of x, the input
%   parameters should correspond to this value of x.
%   - ks corresponds to: [P, M] = ks * [\epsilon, \phi]
n = length(E);
ks = zeros(2, 2);
for i = 1:n
    l = l_vec(YZ(i, :));
    ks = ks + E(i) * A(i) * (l' * l);
end
end
