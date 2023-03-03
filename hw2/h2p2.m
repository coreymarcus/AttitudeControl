%% Setup
clear
close all
clc

addpath("..\util\")

%% Problem Initalization

% Gravity
mu = 398600; %km3/s2

% Orbital elements
a = 6371 + 400;
e = 0;
i = 0;
Ohm = 0;
w = 0;
theta = 270;
[r_inertial, v_inertial] = OE2State(a, e, i, Ohm, w, theta);

% Rotation formulations
J = blkdiag(90, 70, 60); % Inertia matrix
w_b_0 = (pi/180)*[0.0001 0.0001 0]'; % Initial rotation rate, rad/sec
q_inertial2body_0 = [0 0 0 1]'; % Initial attitude quaternion

% Final time
Tf = 60*60;
% Tf = 10000;

% Find orbit rate and assign to omega
n = sqrt(mu/a^3);
w_b_0(3) = -n;

% Nominal rotation
w_b_nominal = [0 0 -n]';

% Angular momentum of the momentum wheel
h0 = 70;

h01 = 0.5*(2*J(1,1)*n + J(2,2)*n + sqrt((-2*J(1,1)*n - J(2,2)*n)^2 - 4*(J(1,1)^2*n^2 + J(1,1)*J(2,2)*n^2)));
h02 = 0.5*(2*J(1,1)*n + J(2,2)*n - sqrt((-2*J(1,1)*n - J(2,2)*n)^2 - 4*(J(1,1)^2*n^2 + J(1,1)*J(2,2)*n^2)));

%% Roll Controller Design

% Parameters
angle_tol = 0.5*pi/180;
safety_fac = 2;

% Maximum roll torque
tau3_max = abs(3*n^2*(J(2,2) - J(1,1)) * (0.5));

% Find natural frequency
omega_n_roll = sqrt(safety_fac*tau3_max/(J(3,3)*angle_tol));
xi_roll = 1;

% Use natural frequency to find proportional gain for roll
% kp_roll = omega_n^2 + 3*n^2*(J(1,1) - J(2,2))/J(3,3)
kp_roll = (omega_n_roll^2 - 3*n^2*(J(1,1) - J(2,2))/J(3,3))*J(3,3)

% Get derivative gain
kd_roll = 2*xi_roll*omega_n_roll*J(3,3)

%% Yaw controller design

% Parameters
angle_tol = 0.5*pi/180;
safety_fac = 2;

% Maximum yaw torque
tau1_max = abs(3*n^2*(J(2,2) - J(3,3)) * (0.5));

% Natural frequency
omega_n_yaw = sqrt(safety_fac*tau1_max/(J(1,1)*angle_tol));
xi_yaw = 1; % For critical damping

% Gains
kp_yaw = J(1,1)*omega_n_yaw^2;
kd_yaw = J(1,1)*2*xi_yaw*omega_n_yaw;

%% Pitch controller design

% Parameters
angle_tol = 0.5*pi/180;
safety_fac = 2;

% Maximum yaw torque
tau2_max = abs(3*n^2*(J(1,1) - J(3,3)) * (0.5));

% Natural frequency
omega_n_pitch = sqrt(safety_fac*tau2_max/(J(2,2)*angle_tol));
xi_yaw = 1; % For critical damping

% Gains
kp_pitch = J(2,2)*omega_n_pitch^2;
kd_pitch = J(2,2)*2*xi_yaw*omega_n_pitch;

%% Main

% Convert OE
[r_inertial_0, v_inertial_0] = OE2State(a, e, i, Ohm, w, theta);

% Run simulation
out = sim("simulink/h2p2_sim");

%% Extract information
r_inertial_hist = squeeze(out.pos);
q_inertial2body_hist = squeeze(out.quat);
w_body_hist = squeeze(out.w);
tout = out.tout;
tau_cont = squeeze(out.tau_cont)';
dtheta = squeeze(out.dtheta);

%% Plotting

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(tout, w_body_hist(ii,:))
    xlabel('Time [s]')
    ylabel('Body Rate [rad/sec]')
end

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(tout, tau_cont(ii,:))
    xlabel('Time [s]')
    ylabel('Control Torque [Nm]')
end

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(tout, dtheta(ii,:))
    xlabel('Time [s]')
    ylabel('Angle Error [rad]')
end