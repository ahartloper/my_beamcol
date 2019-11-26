function [s, Et] = epp_mat(E, sy, e)
% Returns the stress and tangent modulus for Elastic Perf Plastic material.
% @input double E: Elastic modulus
% @input double sy: Yield stress
% @input double e: Strain
% @returns matrix: (double, 2x1) Stress and tangent modulus.
ey = sy / E;
if abs(e) <= ey
    s = E * e;
    Et = E;
else
    s = sign(e) * sy;
    Et = 0.;
end
end
