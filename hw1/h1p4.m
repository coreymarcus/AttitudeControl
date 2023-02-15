%% Setup

clear
close all
clc

addpath("..\util\")

%% Problem Initialization

% Earth radius
Re = 6371;

% Orbital elements
a = Re + 400;
e = 0;
i = 0;
Ohm = 0;
w = 0;
theta = 270;

% Convert OE
[r_inertial, v_inertial] = OE2State(a, e, i, Ohm, w, theta);
mu = 398600; %km3/s2
orbit_rate = sqrt(mu/norm(r_inertial)^3); % Orbit rate, AKA mean motion

% Rotational stuff
J = blkdiag(90,70,60);
w_b_0 = [pi*0.0001/180 pi*0.0001/180 orbit_rate]'; % Initial rotation rate, rad/sec

% Body axes in inertial frame
b_x_in_i = v_inertial/norm(v_inertial);
b_y_in_i = -r_inertial/norm(r_inertial);
b_z_in_i = [0 0 1]';

%% Main

% Setup initial transformation matrix
T_body2inertial0 = [b_x_in_i';
    b_y_in_i';
    b_z_in_i'];
q_body2inertial_0 = DCM2Quat(T_body2inertial0);
q_inertial2body_0 = QuatInv(q_body2inertial_0);

% Propagation timestep
dt = 100;
P = 2*pi/orbit_rate;
Tf = 5*P;

% Time vector
t = 0:dt:Tf;
Nt = length(t);

% Initialize state history vector
state_hist = zeros(13,Nt);
state_hist(1:3,1) = r_inertial;
state_hist(4:6,1) = v_inertial;
state_hist(7:10,1) = q_inertial2body_0;
state_hist(11:13,1) = w_b_0;

state_histGG = state_hist;

% Propagate
for ii = 2:Nt
    state_f = PropagateTwoBody(state_hist(:,ii - 1), dt, J, [0 0 0]',"TwoBody");
    state_hist(:,ii) = state_f;

    state_f_GG = PropagateTwoBody(state_histGG(:,ii - 1), dt, J, [0 0 0]',"TwoBodyAndGG");
    state_histGG(:,ii) = state_f_GG;
end

%% Plotting
figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(t,state_hist(10+ii,:))
    plot(t,state_histGG(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
    legend("No GG", "With GG")
end
