%% Setup
clear
close all
clc

addpath("..\util\")

%% Problem Initalization

% Rotation formulations
J = blkdiag(250, 300, 500); % Inertia matrix
w_b_0 = pi*[1 2 3]'/180; % Initial rotation rate, rad/sec
q_inertial2body_0 = [0 0 0 1]'; % Initial attitude quaternion

% Orbital elements
a = 6371 + 500;
e = 0;
i = 45;
Ohm = 0;
w = 0;
theta = 0;

% Final time
Tf = 6000;

% Controller gain
K_control = 1;

%% Main

% Convert OE
[r_inertial, v_inertial] = OE2State(a, e, i, Ohm, w, theta);

% Propagation timestep
dt = 10;

% Time vector
t = 0:dt:Tf;
Nt = length(t);

% Initialize state history vector
state_hist = zeros(13,Nt);
state_hist(1:3,1) = r_inertial;
state_hist(4:6,1) = v_inertial;
state_hist(7:10,1) = q_inertial2body_0;
state_hist(11:13,1) = w_b_0;

% Propagate
for ii = 2:Nt
    % Find control torque
    T = BdotController(K_control,state_hist(1:3,ii - 1),state_hist(7:10,ii-1), state_hist(11:13,ii - 1));

    state_f = PropagateTwoBody(state_hist(:,ii - 1), dt, J, T);
    state_hist(:,ii) = state_f;
end

%% Plotting
close all

figure
plot3(state_hist(1,:),state_hist(2,:),state_hist(3,:))
axis(1.5*[-a a -a a -a a])
xlabel("x [km]")
ylabel("y [km]")
zlabel("z [km]")
title("Vehicle Position in Inertial Space")

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(t,state_hist(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
end