clear
close all
clc

addpath("..\util\")

%% Problem Initalization

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
Phat0 = blkdiag(1E-3*eye(3),1E-3*eye(3),1E-3*eye(3));

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
    bias_hist_true(:,ii) = bias_hist_true(:,ii-1) + nu_bias(:,ii-1);

    % Generate measurements
    gyro_meas(:,ii) = bias_hist_true(:,ii) + nu_gyro(:,ii);
    B_meas(:,ii) = DipoleMagneticField(state_f(1:3)') + nu_mag(:,ii);
end


%% Plotting

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(t,bias_hist_true(ii,:))
    title("True Bias")
    xlabel("Time [sec]")
    ylabel("Bias [rad/sec]")
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