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


%% Main

% Find orbit rate and assign to omega
n = sqrt(mu/a^3);
w_b_0(3) = n;

% Convert OE
[r_inertial_0, v_inertial_0] = OE2State(a, e, i, Ohm, w, theta);

% Run simulation
out = sim("h2p2_sim");

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