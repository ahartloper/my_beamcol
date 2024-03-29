%% Assignment 5 Solution
%
% This document contains the solution to Assignment 5 of CIVE449 at EPFL. The
% MATLAB source files are available upon request.
%
%% 200x200 constant cross-section and 3,5pt Gauss-Lobatto integration
%

clear all;
close all;
clc;

% Member properties
a0 = 200.;
b0 = 200.;
aL = 200.;
A_0 = a0 * b0;
A_L = aL * b0;
ee = 200000.;
sy0 = 355.;
L = 2000.;
n_int_pts_5 = 5;
n_int_pts_3 = 3;

% Define the fiber section and properties along the length
% Use nf fibers along y, 1 along z
nf = 8;
offset = b0 / nf;
YZ = zeros(nf, 2);
YZ(:, 1) = -b0 / 2 + offset / 2 + offset * (0:(nf-1));

% Area per fiber at integration points
Af_0 = A_0 / nf;
Af_L = A_L / nf;
Af_L2 = (Af_0 + Af_L) / 2.;

% Each column in E and A corresponds to an integration point
% 3pt - 1: r = 0, 2: r = -1, 3: r = 1
% 5pt - 1: r = 0, 2: r = -1, 3: r = 1, 4: r = -sqrt(3./7.), 5: r = sqrt(3./7.)
A5 = Af_L * ones(nf, n_int_pts_5);
A3 = Af_L * ones(nf, n_int_pts_3);

% Constitutive law
c = @(e) epp_mat(ee, sy0, e);
% c = @(e) lin_hard_mat(1., 1., 0.1, e);

% Boundary conditions
fixed_dof = [1, 2, 3];

% Order: A, YZ, L, c, fixed_dof, n_int_pts
beam_model_5 = {A5, YZ, L, c, fixed_dof, n_int_pts_5};
beam_model_3 = {A3, YZ, L, c, fixed_dof, n_int_pts_3};


%% Validate stiffness matrix

E = ee * ones(nf, n_int_pts_5);
S = zeros(nf, n_int_pts_5);
[Ke, Q] = disp_fiber_GL(A5, E, S, YZ, L, n_int_pts_5);

% Validation matrix, see McGuire, Gallagher, and Ziemian pg. 73 (4.32)
I = b0 * a0^3 / 12.;
ke_test = ee * [A_0/L, 0, 0, -A_0/L, 0, 0;
    0, 12*I/L^3, 6*I/L^2, 0, -12*I/L^3, 6*I/L^2;
    0, 6*I/L^2, 4*I/L, 0, -6*I/L^2, 2*I/L;
    -A_0/L, 0, 0, A_0/L, 0, 0;
    0, -12*I/L^3, -6*I/L^2, 0, 12*I/L^3, -6*I/L^2;
    0, 6*I/L^2, 2*I/L, 0, -6*I/L^2, 4*I/L];

error = norm(Ke - ke_test) / norm(ke_test)
% 8 fibers are required for the error to be less than 2%.

%% Validate pure tension

S = ones(nf, n_int_pts_5);
[Ke, Q] = disp_fiber_GL(A5, E, S, YZ, L, n_int_pts_5);
Q_test = [-a0*b0; 0; 0; a0*b0; 0; 0];
error = norm(Q - Q_test) / norm(Q_test)


%% Set the initial conditions and run load increments

V_end = 450.e3;
q0 = [0; 0; 0; 0; 0; 0];
Qfinal = [0; 0; 0; 0; V_end; 0.];
n_increments = 30;
q_incr_5 = load_control(q0, Qfinal, beam_model_5, n_increments);
q_incr_3 = load_control(q0, Qfinal, beam_model_3, n_increments);

%% Validate elastic stiffness

K_theory = 3. * ee * I / L^3;
K_model = q_incr_5{1}(5, 2) / q_incr_5{1}(5, 1);
error = (K_model - K_theory) / K_theory
% Note that the error in the stiffness is similar to the difference in the matrix norms.

%% Plot member response
plot_dof = 5;
figure
hold on;
for i = 1:n_increments
    if i == 1
        plot(q_incr_5{i}(plot_dof, 1), q_incr_5{i}(plot_dof, 2), 'ko', 'displayname', '5pt')
        plot(q_incr_3{i}(plot_dof, 1), q_incr_3{i}(plot_dof, 2), 'ro', 'displayname', '3pt')
    else
        plot(q_incr_5{i}(plot_dof, 1), q_incr_5{i}(plot_dof, 2), 'ko', 'handlevisibility', 'off')
        plot(q_incr_3{i}(plot_dof, 1), q_incr_3{i}(plot_dof, 2), 'ro', 'handlevisibility', 'off')
    end
end
legend('location', 'best')
xlabel(['q_',num2str(plot_dof)]);
ylabel(['Q_',num2str(plot_dof)]);

%% Plot the stresses
i_incr = n_increments;

% 5pt integration
r = [0.; -1.; 1.; -sqrt(3./7.); sqrt(3./7.)];
q = q_incr_5{i_incr}(:, 1);
figure
hold on;
for i = 1:length(r)
    [E, S] = assemble_ES(r(i), q, YZ, L, c);
    plot(S, YZ(:, 1), 'o', 'displayname', num2str(r(i)))
end
xlabel('Stress')
ylabel('Fiber y loc')
legend('location', 'best')

% 3pt integration
r = [0.; -1.; 1.];
q = q_incr_3{i_incr}(:, 1);
figure
hold on;
for i = 1:length(r)
    [E, S] = assemble_ES(r(i), q, YZ, L, c);
    plot(S, YZ(:, 1), 'o', 'displayname', num2str(r(i)))
end
xlabel('Stress')
ylabel('Fiber y loc')
legend('location', 'best')

%% Determine the point of first yield for the 5pt model
figure
hold on;
for i = 1:n_increments
    plot(i, q_incr_5{i}(plot_dof, 2), 'ko', 'handlevisibility', 'off')
end
xlabel('Iteration')
ylabel('Q_5')

% With 8 fiber elements the first yield is between increments 17 and 18 with the
% applied loading scheme.
i_incr = 17;
% 5pt integration
r = [0.; -1.; 1.; -sqrt(3./7.); sqrt(3./7.)];
q = q_incr_5{i_incr}(:, 1);
figure
hold on;
for i = 1:length(r)
    [E, S] = assemble_ES(r(i), q, YZ, L, c);
    plot(S, YZ(:, 1), 'o', 'displayname', num2str(r(i)))
end
xlabel('Stress')
ylabel('Fiber y loc')
legend('location', 'best')

yield_incr = i_incr;
Sx = a0 * b0^2 / 6.;
My = Sx * sy0
Vy = My / L
Qy = q_incr_5{yield_incr}(plot_dof, 2)

error = (Qy - Vy) / Vy

% The difference between the two values is relatively large because there are
% only 8 fibers used for the cross-section. When 50 fibers are used, the values
% are within around 1% of Vy at increment 16.

%% Varying cross-section

% 5pt - 1: r = 0, 2: r = -1, 3: r = 1, 4: r = -sqrt(3./7.), 5: r = sqrt(3./7.)
r = [0.; -1.; 1.; -sqrt(3./7.); sqrt(3./7.)];
width_vary = @(r, l) 1. - 0.5 * (r * l / 2. + l / 2.) / l;
A5_vary = Af_L * ones(nf, n_int_pts_5);
A5_vary(:, 1) = A5_vary(:, 1) * width_vary(r(1), L);
A5_vary(:, 2) = A5_vary(:, 2) * width_vary(r(2), L);
A5_vary(:, 3) = A5_vary(:, 3) * width_vary(r(3), L);
A5_vary(:, 4) = A5_vary(:, 4) * width_vary(r(4), L);
A5_vary(:, 5) = A5_vary(:, 5) * width_vary(r(5), L);

% Order: A, YZ, L, c, fixed_dof, n_int_pts
beam_model_5_vary = {A5_vary, YZ, L, c, fixed_dof, n_int_pts_5};
q_incr_5_vary = load_control(q0, Qfinal, beam_model_5_vary, n_increments);

plot_dof = 5;
figure
hold on;
for i = 1:n_increments
    if i == 1
        plot(q_incr_5{i}(plot_dof, 1), q_incr_5{i}(plot_dof, 2), 'ko', 'displayname', '5pt')
        plot(q_incr_5_vary{i}(plot_dof, 1), q_incr_5_vary{i}(plot_dof, 2), 'bo', 'displayname', '5pt-vary')
    else
        plot(q_incr_5{i}(plot_dof, 1), q_incr_5{i}(plot_dof, 2), 'ko', 'handlevisibility', 'off')
        plot(q_incr_5_vary{i}(plot_dof, 1), q_incr_5_vary{i}(plot_dof, 2), 'bo', 'handlevisibility', 'off')
    end
end
legend('location', 'best')
xlabel(['q_',num2str(plot_dof)]);
ylabel(['Q_',num2str(plot_dof)]);
xlim([0, 250])

% Doesn't converge for the last two increments because the cross-section is
% completely plastic at the support.
i_incr = 28;
q = q_incr_5_vary{i_incr}(:, 1);
figure
hold on;
for i = 1:length(r)
    [E, S] = assemble_ES(r(i), q, YZ, L, c);
    plot(S, YZ(:, 1), 'o', 'displayname', num2str(r(i)))
end
xlabel('Stress')
ylabel('Fiber y loc')
legend('location', 'best')