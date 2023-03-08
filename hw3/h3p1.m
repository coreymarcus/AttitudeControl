%% Setup
clear
close all
clc

addpath("..\util\")

%% Problem Initalization

% Gravity
mu = 398600; %km3/s2

% Orbital elements
a = 6371 + 500000;
e = 0;
i = 0;
Ohm = 0;
w = 0;
theta = 0;
[r_inertial_0, v_inertial_0] = OE2State(a, e, i, Ohm, w, theta);

% Find orbit rate and assign to omega
n = sqrt(mu/a^3);

% Rotation formulations
J = [67946, -83 11129;
    -83 90061 103;
    11129 103 45821];
w_b_0 = [0 0 0]'; % Initial rotation rate, rad/sec
q_inertial2body_0 = [0 0 0 1]'; % Initial attitude quaternion

% Nominal rotation rate
w_b_nom = pi/180*[0 0 1/300]';

% Final simulation time
Tf = 1800;
% Tf = 10000;

% Final manuever time
Tf_man = 900;

% Reaction Wheel Orientation
X_rw_in_body = [0, sqrt(3)/2*cosd(30), sqrt(3)/2*cosd(30), 0, -sqrt(3)/2*cosd(30), -sqrt(3)/2;
    cosd(30), 1/2*cosd(30), -1/2*cosd(30), -cosd(30), -1/2*cosd(30), 1/2*cosd(30);
    sind(30), sind(30), sind(30), sind(30), sind(30), sind(30)];

% Reaction wheel moment
Jw = 0.1295;

%% Roll Controller Design

tau = 0.001;

kp_roll = 2;

% Get derivative gain
kd_roll = kp_roll/tau;

%% Yaw controller design

% Gains
kp_yaw = 2;
kd_yaw = kp_yaw/tau;

%% Pitch controller design

% Gains
kp_pitch = 2;
kd_pitch = kp_pitch/tau;

%% Find final attitude

state0 = zeros(13,1);
state0(1:3) = r_inertial_0;
state0(4:6) = v_inertial_0;
state0(7:10) = q_inertial2body_0;
state0(11:13) = w_b_nom;
statef = PropagateTwoBody(state0,Tf_man,eye(3),[0 0 0]','TwoBody')';

q_inertial2body_f = statef(7:10);

%% Main

out = sim("simulink\h3p1_sim.slx");

%% Extract information
r_inertial_hist = squeeze(out.pos);
q_inertial2body_hist = squeeze(out.quat);
q_inertial2body_ref = squeeze(out.ref_quat);
w_body_hist = squeeze(out.w);
w_body_ref = squeeze(out.ref_rate);
tout = out.tout;
tau_cont = squeeze(out.tau_cont)';
dtheta = squeeze(out.dtheta);
rw_speed = squeeze(out.rw_speed)';

%% Plotting

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(tout, w_body_hist(ii,:))
    plot(tout, w_body_ref(ii,:))
    xlabel('Time [s]')
    ylabel('Body Rate [rad/sec]')
    legend("Actual","Reference")
end

figure
for ii = 1:6
    subplot(6,1,ii)
    hold on
    plot(tout, rw_speed(ii,:))
    xlabel('Time [s]')
    ylabel('Reaction Wheel Speed [rad/sec]')
end

figure
for ii = 1:4
    subplot(4,1,ii)
    hold on
    plot(tout, q_inertial2body_hist(ii,:))
    plot(tout, q_inertial2body_ref(ii,:))
    xlabel('Time [s]')
    ylabel('Quat Element')
    legend("Actual","Reference")
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