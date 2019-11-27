function [q, Qfinal] = nl_analysis_fiber_GL(q0, Qext, fixed_dof, A, YZ, L, c, npts, incr_num)
% Runs Newton iterations for an fiber beam model increment.
% @input matrix q0: (double, 6x1) Initial nodal displacement.
% @input matrix Qext: (double, 6x1) External load vector.
% @input matrix fixed_dof: (int, 1xr) Fixed degrees of freedom, 1 <= r_i <= 6,
%   no values repeated.
% @input matrix A: (double, nxm) Area of each fiber at the integration points.
% @input matrix YZ: (double, nx2) y,z-coordinates of each fiber.
% @input double L: Element length.
% @input function c: Constitutive law, stress = c(strain).
% @input int npts: Number of Gauss-Lobatto integration points.
% @input incr_num: Current increment.
%
% Notes:
%   - Uses fiber element with Gauss-Lobatto integration
%   - Assumes constant YZ for all the integration points
%   - Assumes that A does not change throughout the analysis
%   - Newton tolerance is set on the norm of the relative force residual to 1e-8
%   - Maximum of 16 iterations are run

% Initalization
tolerance = 1.e-8;
max_it = 16;
free_dof = 1:6;
free_dof(fixed_dof) = [];
if npts == 3
    r = [0.; -1.; 1.];
elseif npts == 5
    r = [0.; -1.; 1.; -sqrt(3./7.); sqrt(3./7.)];
end

% Pre iteration 
q = q0;
resid = zeros(6, 1);
[E, S] = assemble_ES(r, q, YZ, L, c);
[Ke, Qint] = disp_fiber_GL(A, E, S, YZ, L, npts);
resid(free_dof) = Qext(free_dof) - Qint(free_dof);

% Newton iterations for the load step
for i = 1:max_it
    kff = Ke(free_dof, free_dof);
    Qff = resid(free_dof);
    qff = kff \ Qff;
    q(free_dof) = q(free_dof) + qff;
    [E, S] = assemble_ES(r, q, YZ, L, c);
    [Ke, Qint] = disp_fiber_GL(A, E, S, YZ, L, npts);
    resid(free_dof) = Qext(free_dof) - Qint(free_dof);
    conv_criteria = norm(resid) / norm(Qext);
    if conv_criteria < tolerance
        break;
    end
end
fprintf('%d\t\t%d\t\t\t%0.3e\n', incr_num, i, conv_criteria);
if i == max_it
    fprintf('Increment %d did not converge!\n', incr_num);
end

Qfinal = Qext;
Qfinal(fixed_dof) = Ke(fixed_dof, free_dof) \ q(free_dof);

end
