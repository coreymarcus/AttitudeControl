function [B_i] = DipoleMagneticField(r_inertial)
%DipoleMagneticField finds the magentic field strength (B_i) in the
%inertial frame as a function of and inertial, earth centered position,
%r_inertial [km]. The function is hard coded for Earth.

% Constants
B0 = 3.12E-5; % Mean magnetic field at equater, [Tesla]
Re = 6371; % Earth radius [km]

% Find latitude of spacecraft
r_xy = norm(r_inertial(1:2));
lambda = atan2(r_inertial(3),r_xy);

% Find B in NED
B_n = B0*(Re/norm(r_inertial))^3*[cos(lambda); 0; 2*sin(lambda)];

% Transformation to inertial space
T_inertial2NED = Inertial2NED(r_inertial);

% Output
B_i = T_inertial2NED'*B_n;

end