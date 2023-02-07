function dxdt = AttitudeEvolution(~,q_omega,J,T)
%AttitudeEvolution finds the instantaneous rate of change of q_omega =
%[quat_inertial2body; omega_bodyWRTinerital_inbody] as a function of
%inertia matrix J and input torque T

% Extract locals
q = q_omega(1:4);
omega = q_omega(5:7);

% Change in omega
omega_dot = J\(T - cross(omega, J*omega));
q_dot = 0.5*QuatProduct([omega; 0],q);

% Output
dxdt = [q_dot; omega_dot];

end