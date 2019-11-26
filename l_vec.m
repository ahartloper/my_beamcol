function l = l_vec(yz)
% Returns the geometrical vector for a single fiber.
% @param matrix yz: (double, 1x2) y,z coordinates of the fiber
% @returns matrix: (double, 1x2) Geometrical vector
%
l = [1., yz(1)];
end