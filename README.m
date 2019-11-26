%% Instructions
% 
% This MATLAB library provides a displacement beam-column element that is
% integrated using Gauss-Lobatto integration (3 and 5pt options are available).
% A load-controlled analysis is provided, only a single element can be
% analyzed. The nonlinear problem is solved using Newton's method. Varying
% area can be considered, but the location of the fiber sections are fixed
% across all the integration points. Nonlinear material can be considered.
% Changes in the geometry are not considered.
%  
% 
%% Define the beam model
%
% The beam model is defined by: 
%
% 
% * A: Area of each fiber at the integration points
% * YZ: (y,z) Coordinates for all of the fibers
% * L: Length of the element
% * c: Constitutive model for the fiber elements
% * fixed_dof: Fixed DOF in the analysis (boundary conditions)
% * npts: Number of integration points (Gauss-Lobatto)
%
% Then set:
%
%
%   model = {A, YZ, L, c, fixed_dof, npts}
% 
%
%% Set the initial displacement, external load, and increments
% 
% The initial displacement is typically taken as a vector of zeros (no
% deformation).
% The external load is simply defined by the load at the end of the step. The
% load control will divide the external load evenly by the number of increments.
% 
% An example is:
%
%
%   q0 = [0; 0; 0; 0; 0; 0];
%   Qfinal = [0; 0; 0; 0; 0; 1];
%   n_increments = 10;
% 
%



%% Run the load control
% 
% The analysis is simply run using
%
%
%   q_incr = load_control(q0, Qfinal, model, n_increments);
%
%
% Where all the arguments to this function have been defined. A cell array is
% output, where the nodal displacements and nodal forces are provided at each
% increment. Note that the reaction forces are not calculated.
