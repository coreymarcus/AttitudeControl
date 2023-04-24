function [qf, wf] = AttitudePropagate(q0, w0, J, dt, method, T)
%AttitudePropagate propagates a quaternion for time dt

odefun = @(t, state) AttitudeEvolution(t,state,J,T);

% ode options
switch method
    case "ode45"
        odeoptions = odeset("RelTol",1E-8,"AbsTol",1E-8);

        % Propagate
        [~, stateout] = ode45(odefun,[0 dt],[q0; w0], odeoptions);
        N_state = size(stateout,1);

    case "RK4"
        N_state = 100;
        dt_RK4 = dt/(N_state-1);
        [stateout, ~] = RK4(0, dt_RK4, N_state, [q0; w0], odefun);

    otherwise
        error("Invalid Propagation Method")
end

% Final state
state_f = stateout(N_state,1:7)';
qf = state_f(1:4);
wf = state_f(5:7);

% Renormalize final state
qf = qf/norm(qf);

end