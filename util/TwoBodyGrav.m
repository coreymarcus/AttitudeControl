function dxdt = TwoBodyGrav(~,rv_inertial)
% TwoBodyGrav finds the instantaneous rate of change for a spacecraft at
% rv_inertial = [r_inertial; v_inertial] [km] and [km/sec]. Hard coded for
% Earth

mu = 398600; %km3/s2

%extract state information
r_x = rv_inertial(1);
r_y = rv_inertial(2);
r_z = rv_inertial(3);
v_x = rv_inertial(4);
v_y = rv_inertial(5);
v_z = rv_inertial(6);

%calculate distance from center of earth
nr = norm([r_x r_y r_z]);

%calculate accelerations
a_x = -mu*r_x/nr^3;
a_y = -mu*r_y/nr^3;
a_z = -mu*r_z/nr^3;

%generate output vector
dxdt = [v_x v_y v_z a_x a_y a_z]';

end