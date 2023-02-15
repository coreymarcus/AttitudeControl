%% Setup
clear
close all
clc

addpath("..\util\")

%% Problem Initalization

% Rotation formulations
J = blkdiag(60, 90, 60); % Inertia matrix
w_b_0 = [0 0.5 0]'; % Initial rotation rate, rad/sec
q_inertial2body_0 = [0 0 0 1]'; % Initial attitude quaternion

% Orbital elements
a = 6371 + 500;
e = 0;
i = 0;
Ohm = 0;
w = 0;
theta = 0;

% Final time
Tf = 60*60;

% Gravity
mu = 398600; %km3/s2

% Constant input torque
T = [0.001 0 0]';

%% Main

% Convert OE
[r_inertial_0, v_inertial_0] = OE2State(a, e, i, Ohm, w, theta);

% Run simulation
out = sim("h2p1_sim");

%% Extract information
r_inertial_hist = squeeze(out.pos);
q_inertial2body_hist = squeeze(out.quat);
w_body_hist = squeeze(out.w);
tout = out.tout;

%% Plotting

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(tout, w_body_hist(ii,:))
    xlabel('Time [s]')
    ylabel('Body Rate [rad/sec]')
end
