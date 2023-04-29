clear
close all
clc

addpath("..\util\","..\project")

%% Problem Initalization

rng(4)

% Rotation formulations
J = blkdiag(1, 1, 1); % Inertia matrix
w_b_0 = [0 0 0]'; % Initial rotation rate, rad/sec
q_inertial2body_0 = [0 0 0 1]'; % Initial attitude quaternion

% Orbital elements
a = 6371 + 400;
e = 0;
i = 45;
Ohm = 0;
w = 0;
theta = 0;
mu = 398600; %km3/s2

% Final time
Tp = 2*pi*sqrt(a^3/mu);
Tf = 8*Tp;

% Propagation timestep
dt = 1;

% Controller gain
K_control = 1;

% Variance for magnetometer noise
Sigma_mag = (2.5E-7)^2*eye(3);

% Noise characteristics for IMU
gyro_ARW = 0.15; % deg/sqrt(hour)
bias_inst = 0.3; % deg/hour

% Initial estimate uncertainty
Phat0 = blkdiag(1E-3*eye(3),1E-10*eye(3));

%% Main

% Convert IMU characteristics to something usable
Sigma_nu = ((gyro_ARW*(pi/180))^2)*(1/3600)*dt*eye(3); %(([deg/sqrt(hour)]*[rad/deg])^2)*[hour/sec]*[sec]
Sigma_b = (bias_inst*(pi/180)*(1/3600)*dt)^2*eye(3); %([deg/hour]*[rad/deg]*[hour/sec]*[sec])^2

% Convert OE
[r_inertial, v_inertial] = OE2State(a, e, i, Ohm, w, theta);

% Time vector
t = 0:dt:Tf;
Nt = length(t);

% Initialize state history vector
state_hist_true = zeros(13,Nt);
state_hist_true(1:3,1) = r_inertial;
state_hist_true(4:6,1) = v_inertial;
state_hist_true(7:10,1) = q_inertial2body_0;
state_hist_true(11:13,1) = w_b_0;

% Estimate histories
bias_est = zeros(3,Nt);
q_est = zeros(4,Nt);
Phat = zeros(6,6,Nt);
q_est_err = zeros(3,Nt);

% Initial estimates
init_est = mvnrnd(zeros(1,6),Phat0)';
bias_est(:,1) = init_est(4:6);
dq0 = [0.5*init_est(1:3); 1];
dq0 = dq0/norm(dq0);
q_est(:,1) = QuatProduct(dq0,q_inertial2body_0);
Phat(:,:,1) = Phat0;
q_est_err(:,1) = 2*q_est(1:3,1);

% Measurement histories
gyro_meas = zeros(3,Nt);
B_meas = zeros(3,Nt);

% Initialize bias
bias_hist_true = zeros(3,Nt);

% Sample noise
nu_bias = mvnrnd(zeros(1,3),Sigma_b,Nt)';
nu_gyro = mvnrnd(zeros(1,3),Sigma_nu,Nt)';
nu_mag = mvnrnd(zeros(1,3),Sigma_mag,Nt)';

% Propagate
for ii = 2:Nt

    % Propagate truth
    state_f = PropagateTwoBody(state_hist_true(:,ii - 1), dt, J, [0 0 0]',"TwoBodyNoAttitude");
    state_hist_true(:,ii) = state_f;
    bias_hist_true(:,ii) = bias_hist_true(:,ii-1) + 0*nu_bias(:,ii-1);

    % Generate measurements
    gyro_meas(:,ii) = bias_hist_true(:,ii) + nu_gyro(:,ii);
    B_inertial = DipoleMagneticField(state_f(1:3)');
    B_body = NativeQuatTransform(q_inertial2body_0,B_inertial);
    B_meas(:,ii) = B_body + nu_mag(:,ii);

    % Estimate
    [bias_est(:,ii), Phat(:,:,ii), q_est(:,ii)] = MEKF(...
        bias_est(:,ii-1),...
        Phat(:,:,ii-1),...
        dt,...
        B_meas(:,ii),...
        Sigma_mag,...
        q_est(:,ii-1),...
        Sigma_nu,...
        Sigma_b,...
        gyro_meas(:,ii),...
        J,...
        state_hist_true(1:3,ii));
%     [Phat(1:3,1:3,ii), q_est(:,ii)] = MEKF_no_bias(...
%         Phat(1:3,1:3,ii-1),...
%         dt,...
%         B_meas(:,ii),...
%         Sigma_mag,...
%         q_est(:,ii-1),...
%         Sigma_nu,...
%         gyro_meas(:,ii),...
%         J,...
%         state_hist_true(1:3,ii));

    q_est_err(:,ii) = 2*q_est(1:3,ii);
end


%% Plotting

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(t,bias_hist_true(ii,:))
    plot(t,bias_est(ii,:))
    title("Bias")
    xlabel("Time [sec]")
    ylabel("Bias [rad/sec]")
    legend("true",'est')
end

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(t,gyro_meas(ii,:))
    title("Gyro Meas")
    xlabel("Time [sec]")
    ylabel("Gyro Meas [rad/sec]")
end

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(t,B_meas(ii,:))
    title("Mag Field Meas")
    xlabel("Time [sec]")
    ylabel("Mag Field [magnets/sec]")
end

figure
for ii = 1:4
    subplot(4,1,ii)
    plot(t,q_est(ii,:))
    title("Estimated Attitude")
    xlabel("Time [sec]")
    ylabel("Quat Element")
end

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(t,q_est_err(ii,:))
    plot(t,3*sqrt(squeeze(Phat(ii,ii,:))))
    plot(t,-3*sqrt(squeeze(Phat(ii,ii,:))))
    title("Quat Est Error")
    xlabel("Time [sec]")
    ylabel("Quat Element")
end