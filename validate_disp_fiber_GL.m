clear all;
close all;
clc;

TEST_TOL = 1.e-8;

%% 2x2 constant cross-section and 3pt Gauss-Lobatto

% Member properties
a0 = 2;
b0 = 2;
aL = 2.;
A_0 = a0 * b0;
A_L = aL * b0;
ee = 1.;
L = 2.;
n_int_pts = 3;

% Define the fiber section and properties along the length
% Use nf fibers along y, 1 along z
nf = 4;
offset = b0 / nf;
YZ = zeros(nf, 2);
YZ(:, 1) = -b0 / 2 + offset / 2 + offset * (0:(nf-1));

% Area per fiber at integration points
Af_0 = A_0 / nf;
Af_L = A_L / nf;
Af_L2 = (Af_0 + Af_L) / 2.;

% Each column in E and A corresponds to an integration point
% 3pt - 1: r = 0, 2: r = -1, 3: r = 1
% 5pt - 1: r = 0, 2: r = -1, 3: r = 1, 4: r = -sqrt(21), 5: r = sqrt(21)
E = ee * ones(nf, n_int_pts);
A = Af_L * ones(nf, n_int_pts);


%% Validate stiffness matrix

% Use Gauss-Lobatto integration
S = zeros(nf, n_int_pts);
[Ke, Q] = disp_fiber_GL(A, E, S, YZ, L, n_int_pts)

% Validation matrix, see McGuire, Gallagher, and Ziemian pg. 73 (4.32)
I = b0 * a0^3 / 12.;
ke_test = [A_0/L, 0, 0, -A_0/L, 0, 0;
    0, 12*I/L^3, 6*I/L^2, 0, -12*I/L^3, 6*I/L^2;
    0, 6*I/L^2, 4*I/L, 0, -6*I/L^2, 2*I/L;
    -A_0/L, 0, 0, A_0/L, 0, 0;
    0, -12*I/L^3, -6*I/L^2, 0, 12*I/L^3, -6*I/L^2;
    0, 6*I/L^2, 2*I/L, 0, -6*I/L^2, 4*I/L];

error = norm(Ke - ke_test) / norm(ke_test)

%% Validate pure tension

S = ones(nf, n_int_pts);

% Use Gauss-Lobatto integration
[Ke, Q] = disp_fiber_GL(A, E, S, YZ, L, n_int_pts);
Q_test = [-4.; 0; 0; 4.; 0; 0];
error = norm(Q - Q_test) / norm(Q_test);
assert(error < TEST_TOL);

%% Validate pure bending

S = zeros(nf, n_int_pts);
S([1, 4], 1) = [-1, 1];
S([1, 4], 2) = [-1., 1.];
S([1, 4], 3) = [-1., 1.];

% Use Gauss-Lobatto integration
[Ke, Q] = disp_fiber_GL(A, E, S, YZ, L, n_int_pts);
Q_test = [0; 0; -1.5; 0; 0; 1.5];
error = norm(Q - Q_test) / norm(Q_test);
assert(error < TEST_TOL);

%% Validate shear + bending

S = zeros(nf, n_int_pts);
S([1, 4], 1) = [-0.5, 0.5];
S([1, 4], 2) = [-1., 1.];
S([1, 4], 3) = [0., 0.];

% Use Gauss-Lobatto integration
[Ke, Q] = disp_fiber_GL(A, E, S, YZ, L, n_int_pts);
Q_test = [0; -0.75; -1.5; 0; 0.75; 0];
error = norm(Q - Q_test) / norm(Q_test);
assert(error < TEST_TOL);


%% Error in stiffness matrix with respect to the number of fibers along height

error = [];
nf = 4:10:200;
for i = 1:length(nf)
    Af_0 = A_0 / nf(i);
    Af_L = A_L / nf(i);
    Af_L2 = (Af_0 + Af_L) / 2.;
    offset = b0 / nf(i);
    YZ = zeros(nf(i), 2);
    E = ee * ones(nf(i), n_int_pts);
    A = Af_L * ones(nf(i), n_int_pts);
    S = ones(nf(i), n_int_pts);
    YZ(:, 1) = -b0 / 2 + offset / 2 + offset * (0:(nf(i)-1));
    [Ke, Q] = disp_fiber_GL(A, E, S, YZ, L, n_int_pts);

    error = [error; norm(Ke - ke_test) / norm(ke_test)];
end

figure
semilogy(nf, error)
xlabel('Num Fibers')
ylabel('Error')
