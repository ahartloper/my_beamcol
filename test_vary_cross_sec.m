clear all;
close all;
clc;

%% 2x2 constant cross-section and Gauss-Lobatto

% Member properties
a0 = 2;
b0 = 2;
aL = 0.5 * 2.;
A_0 = a0 * b0;
A_L = aL * b0;
ee = 1.;
L = 2.;
n_int_pts_3 = 3;

% Define the fiber section and properties along the length
% Use nf fibers along y, 1 along z
nf = 20;
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
A3 = zeros(nf, n_int_pts_3);
A3(:, 1) = Af_L2;
A3(:, 2) = Af_0;
A3(:, 3) = Af_L;

% Constitutive law
c = @(e) epp_mat(1., 1., e);
% c = @(e) lin_hard_mat(1., 1., 0.1, e);

% Boundary conditions
fixed_dof = [1, 2, 3];

% Order: A, YZ, L, c, fixed_dof, n_int_pts
beam_model_3 = {A3, YZ, L, c, fixed_dof, n_int_pts_3};

%% Set the initial conditions and run load increments

q0 = [0; 0; 0; 0; 0; 0];
Qfinal = [0; 0; 0; 0.; 0.; 2 * 4/6];
n_increments = 50;
q_incr_3 = load_control(q0, Qfinal, beam_model_3, n_increments);

%% Plot member response
plot_dof = 6;
figure
% Plot Q vs q
subplot(2, 1, 1)
hold on;
for i = 1:n_increments
    if i == 1
%         plot(q_incr_5{i}(plot_dof, 1), q_incr_5{i}(plot_dof, 2), 'ko', 'displayname', '5pt')
        plot(q_incr_3{i}(plot_dof, 1), q_incr_3{i}(plot_dof, 2), 'ro', 'displayname', '3pt')
    else
%         plot(q_incr_5{i}(plot_dof, 1), q_incr_5{i}(plot_dof, 2), 'ko', 'handlevisibility', 'off')
        plot(q_incr_3{i}(plot_dof, 1), q_incr_3{i}(plot_dof, 2), 'ro', 'handlevisibility', 'off')
    end
end
legend('location', 'best')
xlabel(['q_',num2str(plot_dof)]);
ylabel(['Q_',num2str(plot_dof)]);

% Plot Q vs increments
subplot(2, 1, 2)
hold on;
for i = 1:n_increments
    if i == 1
%         plot(q_incr_5{i}(plot_dof, 1), q_incr_5{i}(plot_dof, 2), 'ko', 'displayname', '5pt')
        plot(i, q_incr_3{i}(plot_dof, 2), 'ro', 'displayname', '3pt')
    else
%         plot(q_incr_5{i}(plot_dof, 1), q_incr_5{i}(plot_dof, 2), 'ko', 'handlevisibility', 'off')
        plot(i, q_incr_3{i}(plot_dof, 2), 'ro', 'handlevisibility', 'off')
    end
end
legend('location', 'best')
xlabel('Increment');
ylabel(['Q_',num2str(plot_dof)]);

%% Plot the stresses
% i_incr = 29;
i_incr = 50;

% 5pt integration
% r = [0.; -1.; 1.; -sqrt(3./7.); sqrt(3./7.)];
% q = q_incr_5{i_incr}(:, 1);
% figure
% hold on;
% for i = 1:length(r)
%     [E, S] = assemble_ES(r(i), q, YZ, L, c);
%     plot(S, YZ(:, 1), 'o', 'displayname', num2str(r(i)))
% end
% xlabel('Stress')
% ylabel('Fiber y loc')
% legend('location', 'best')

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