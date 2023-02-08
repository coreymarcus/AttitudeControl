function T = BdotController(K, r_inertial, q_inertial2body, w_body)
%BdotController outputs the input torque in the body frame as a function of
%the Bdot gain, K, vehicle position, r, in inertial coordinates, and
%vehicle's angular rate, w, in body frame with respect to the inertial
%coordinates, and quaternion from inertial to body

% Find magnetic field
B_inertial = DipoleMagneticField(r_inertial);
B_body = Quat2DCM(q_inertial2body)*B_inertial;

% Torque
T = (B_body*B_body'/norm(B_body)^2 - eye(3))*K*w_body;

end