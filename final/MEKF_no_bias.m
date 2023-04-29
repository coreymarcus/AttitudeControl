function [Phat, qhat] = MEKF_no_bias(Phat0, dt, B_meas, Sigma_B, q0, Q_w, w_meas, J, r_inertial)

% Predicted omega
w_bar = w_meas;

% Propagate Attitude
[q_bar, ~] = AttitudePropagate(q0, w_bar, J, dt, "RK4",[0 0 0]');

% Propagate state estimate
Aeval = PaPaFunc(dt,...
    w_bar(1),...
    w_bar(2),...
    w_bar(3)+1E-10); % Add very small peturbation to prevent divide by zero
Beval = PaPdomegaFunc(dt,...
    w_bar(1),...
    w_bar(2),...
    w_bar(3)+1E-10); % Add very small peturbation to prevent divide by zero
F = Aeval;
G = Beval;
Pbar = F*Phat0*F' + G*Q_w*G';

% Find Jacobian
Ceval = partialBwrtA(0,0,0,...
    q_bar(1),q_bar(2),q_bar(3),q_bar(4), ...
    r_inertial(1),r_inertial(2),r_inertial(3));
H = Ceval;

% Predict measurement
B_inertial_bar = DipoleMagneticField(r_inertial);
B_meas_bar = NativeQuatTransform(q_bar,B_inertial_bar);

% Kalman gain
K = Pbar*H'/(H*Pbar*H' + Sigma_B);
xhat = K*(B_meas - B_meas_bar);

% Update Covariance
Phat = (eye(3) - K*H)*Pbar*(eye(3) - K*H)' + K*Sigma_B*K';

% Apply update to quaternion
dqhat = [0.5*xhat(1:3); 1];
dqhat = dqhat/norm(dqhat);
qhat = QuatProduct(dqhat,q_bar);

end