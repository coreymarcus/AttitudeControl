function [state_f] = PropagateTwoBody(state_0, dt, J, T)
%PropagateTwoBody propagates an state over time dt using two body gravity.
%state = [position_inertial; velocity_inertial; quat_inertial2body;
%omega_bodyWRTinertial_inbody]. J = inertia matrix. T = constant input torque.

% Construct function for ode45
odefun = @(t, state) [TwoBodyGrav(t, state(1:6)); AttitudeEvolution(t,state(7:13),J,T)];

% ode options
odeoptions = odeset("RelTol",1E-8,"AbsTol",1E-8);

% Propagate
[~, stateout] = ode45(odefun,[0 dt],state_0, odeoptions);

% Final state
state_f = stateout(end,:);

% Renormalize final state
state_f(7:10) = state_f(7:10)/norm(state_f(7:10));

end