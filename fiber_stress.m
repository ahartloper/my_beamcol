function [stress, Et] = fiber_stress(q, B, yz, constitutive_law)
% Returns the stress and tangent modulus for a fiber.
% @input matrix q: (double, 6x1) Section deformations.
% @input matrix B: (double, 2x6) Strain-displacement matrix.
% @input matrix yz: (double, 1x2) y,z-coordinates of the fiber.
% @input function constitutive_law: Provides the strain and tangent modulus for
%   a given strain. 
% @returns matrix: (double, 1x2) Stress and tangent modulus for the fiber.
%
% Notes:
%   - constitutive_law needs to be a function of the strain only.
d = B * q;
strain = d(1) + d(2) * yz(1);
[stress, Et] = constitutive_law(strain);
end
