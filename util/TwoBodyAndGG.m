function dxdt = TwoBodyAndGG(t, state, J, T, n)
%TwoBodyAndGG Finds the instantaneous rate of change in the state due to
%two body gravity, input torques, and gravity gradient effects

% Local down direction
q_inertial2body = state(7:10);
down_inertial = -state(1:3)/norm(state(1:3));
down_body = QuatTransform(q_inertial2body,down_inertial);

% Gravity gradient torque
T_GG = 3*n^2*cross(down_body,J*down_body);

% Dynamics
dxdt = [TwoBodyGrav(t, state(1:6)); AttitudeEvolution(t,state(7:13),J,T+T_GG)];
end