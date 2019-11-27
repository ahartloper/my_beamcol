function F = gl_quad(f, npts, L)
% Returns the integral of f(r) using Gauss-Lobatto qudrature.
% @param function f: Integrand as a function of r only
% @param int npts: Number of integration points to use
% @param double L: Original length of the element
% @returns: Integral of f.
%
% Notes:
%   - Only 3 and 5 pt integration supported at the moment.
%   - Integration done in natural coordinate system, -1 < r < 1
%   - f can be made to be a function of r using an anonymous function
if npts == 3
    r = [0.; -1.; 1.];
    weights = [4. / 3.; 1. / 3.; 1. / 3.];
elseif npts == 5
    r = [0.; -1.; 1.; -sqrt(3./7.); sqrt(3./7.)];
    weights = [32./45; 1./10.; 1./10.; 49./90; 49./90.];
end

F = weights(1) * f(r(1));
for i = 2:npts
    F = F + weights(i) * f(r(i));
end
F = L / 2. * F;
end
