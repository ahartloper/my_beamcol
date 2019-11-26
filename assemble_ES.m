function [E, S] = assemble_ES(r, q, YZ, L, c_law)
% Provides the tangent modulus and stress for all fibers at all integration points.
% @input matrix r: (double, mx1) Integration points in natural coordinates
% @input matrix q: (double, 6x1) Nodal deformations
% @input matrix YZ: (double, nx2) y,z-coordinates of each fiber
% @input double L: Element length
% @input function c_law: Constitutive law, stress = f(strain).
% @returns matrix: (double, nxm/nxm) Stress and tangent modulus for all the fibers.
%
% Notes:
%   - constitutive_law needs to be a function of the strain only.
npts = length(r);
nf = length(YZ);
E = zeros(nf, npts);
S = zeros(nf, npts);
% Compute stress for each fiber at each integration point
for j = 1:length(r)
    x = r(j) * L / 2. + L / 2.;
    Bx = b_matrix(x, L);
    for i = 1:nf
        [S(i, j), E(i, j)] = fiber_stress(q, Bx, YZ(i, :), c_law);
    end
end
end
