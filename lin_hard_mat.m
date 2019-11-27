function [s, Et] = lin_hard_mat(E, sy, a, e)
% Returns the stress and tangent modulus for linear hardening material.
% @input double E: Elastic modulus
% @input double sy: Yield stress
% @input double a: Strain hardening ratio
% @input double e: Strain
% @returns matrix: (double, 2x1) Stress and tangent modulus.
ey = sy / E;
if abs(e) <= ey
    s = E * e;
    Et = E;
else
    s = sign(e) * (sy + (abs(e) - ey) * E * a);
    Et = E * a;
end
end
