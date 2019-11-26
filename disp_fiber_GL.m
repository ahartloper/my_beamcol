function [ke, Q] = disp_fiber_GL(A, E, S, YZ, L, npts)
% Returns the element stiffness matrix for a fiber beam-column element.
% @param matrix A: (double, nx3) Area of each fiber at each integration point
% @param matrix E: (double, nx3) Elastic modulus of each fiber at each integration point
% @param matrix S: (double, nx3) Stress of each fiber at each integration point
% @param matrix YZ: (double, nx2) y,z coordinates of each fiber
% @param double L: Length of the element.
% @param int npts: Number of integration points along the length.
% @returns matrix: (double, 6x6) Stiffness matrix for the element.
%
% Notes:
%   - The A, E, S must be provided with respect to the natural cooridinate system
%   - The columns of A, E, S are ordered as:
%       1: r = 0, 2: r = -1, 3: r = 1, ... additional points
%   - This function assumes that YZ is constant along the length
%   - Uses Gauss-Lobatto quadrature to integrate along the length
%   - Only 3 and 5 pt integration is supported at the moment
if npts == 3
    kfun = @(r) k_3pt(r, A, E, YZ, L);
    ffun = @(r) f_3pt(r, A, S, YZ, L);
elseif npts == 5
    kfun = @(r) k_5pt(r, A, E, YZ, L); 
    ffun = @(r) f_5pt(r, A, S, YZ, L);   
end
ke = gl_quad(kfun, npts, L);
Q = gl_quad(ffun, npts, L);
end

%% Helper functions
function k = k_3pt(r, A, E, YZ, L)
% Function to help out with the 3pt integration for varying A and E.
if r == 0.
    ks = assemble_ks(E(:, 1), A(:, 1), YZ);
elseif r == -1.
    ks = assemble_ks(E(:, 2), A(:, 2), YZ);
elseif r == 1.
    ks = assemble_ks(E(:, 3), A(:, 3), YZ);
end
% B is calculated using x = r * L/2 + L/2
x = r * L / 2. + L / 2.;
B = b_matrix(x, L);
k = B' * ks * B;
end

function Q = f_3pt(r, A, sigma, YZ, L)
% Function to help out with the 3pt integration for varying sigma and E.
if r == 0.
    fs = assemble_fs(sigma(:, 1), A(:, 1), YZ);
elseif r == -1.
    fs = assemble_fs(sigma(:, 2), A(:, 2), YZ);
elseif r == 1.
    fs = assemble_fs(sigma(:, 3), A(:, 3), YZ);
end
% B is calculated using x = r * L/2 + L/2
x = r * L / 2. + L / 2.;
B = b_matrix(x, L);
Q = B' * fs;
end


function k = k_5pt(r, A, E, YZ, L)
% Function to help out with the 5pt integration for varying A and E.
if r == 0.
    ks = assemble_ks(E(:, 1), A(:, 1), YZ);
elseif r == -1.
    ks = assemble_ks(E(:, 2), A(:, 2), YZ);
elseif r == 1.
    ks = assemble_ks(E(:, 3), A(:, 3), YZ);
elseif r == -sqrt(3./7.)
    ks = assemble_ks(E(:, 4), A(:, 4), YZ);
elseif r == sqrt(3./7.)
    ks = assemble_ks(E(:, 5), A(:, 5), YZ);
end
% B is calculated using x = r * L/2 + L/2
x = r * L / 2. + L / 2.;
B = b_matrix(x, L);
k = B' * ks * B;
end

function Q = f_5pt(r, A, sigma, YZ, L)
% Function to help out with the 3pt integration for varying sigma and E.
if r == 0.
    fs = assemble_fs(sigma(:, 1), A(:, 1), YZ);
elseif r == -1.
    fs = assemble_fs(sigma(:, 2), A(:, 2), YZ);
elseif r == 1.
    fs = assemble_fs(sigma(:, 3), A(:, 3), YZ);
elseif r == -sqrt(3./7.)
    fs = assemble_fs(sigma(:, 4), A(:, 4), YZ);
elseif r == sqrt(3./7.)
    fs = assemble_fs(sigma(:, 5), A(:, 5), YZ);
end
% B is calculated using x = r * L/2 + L/2
x = r * L / 2. + L / 2.;
B = b_matrix(x, L);
Q = B' * fs;
end
