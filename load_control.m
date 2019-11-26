function q_results = load_control(q0, Qfinal, model, n_increments)
% Runs load control analysis for a fiber beam model.
% @input matrix q0: (double, 6x1) Initial nodal displacement.
% @input matrix Qfinal: (double, 6x1) Total applied external load.
% @input cell model: Properties for the beam model.
% @input int n_incrments: Number of increments to divide Qfinal.
% @returns cell: {q_i, Q_i} Nodal displacements and forces at each increment.
% 
% Notes:
%   - model has the order: A, YZ, L, c, fixed_dof, n_int_pts
A = model{1};
YZ = model{2};
L = model{3};
c = model{4};
fixed_dof = model{5};
npts = model{6};
% Run the analysis
q = q0;
q_results = {};
fprintf('Incr.\tIts.\t\t||Resid.||\n');
for incr_num = 1:n_increments
    Qext = incr_num / n_increments * Qfinal;
    [q, Qf] = nl_analysis_fiber_GL(q, Qext, fixed_dof, A, YZ, L, c, npts, incr_num);
    q_results = [q_results, [q, Qf]];
end
