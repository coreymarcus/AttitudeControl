%% Setup
clear
close all
clc

addpath("..\util\")

%% Problem Initalization

% Rotation formulations
J = blkdiag(100, 75, 50); % Inertia matrix
w_b_0_1 = [0 .1 0]'; % Initial rotation rate, rad/sec
w_b_0_2 = [0.01 .0001 0.0001]'; % Initial rotation rate, rad/sec
w_b_0_3 = [0.0001 .01 0.0001]'; % Initial rotation rate, rad/sec
q_inertial2body_0 = [0 0 0 1]'; % Initial attitude quaternion

% We don't actually need position and velocity but easier to interface with
% the code this way
a = 6371 + 500;
e = 0;
i = 0;
Ohm = 0;
w = 0;
theta = 0;

% Final time
Tf = 2.5*60*60;

%% Main

% Convert OE
[r_inertial, v_inertial] = OE2State(a, e, i, Ohm, w, theta);

% Propagation timestep
dt = 10;

% Time vector
t = 0:dt:Tf;
Nt = length(t);

% Initialize state history vector
state_hist1 = zeros(13,Nt);
state_hist1(1:3,1) = r_inertial;
state_hist1(4:6,1) = v_inertial;
state_hist1(7:10,1) = q_inertial2body_0;
state_hist2 = state_hist1;
state_hist3 = state_hist1;
state_hist1(11:13,1) = w_b_0_1;
state_hist2(11:13,1) = w_b_0_2;
state_hist3(11:13,1) = w_b_0_3;

% Propagate
for ii = 2:Nt
    state_f1 = PropagateTwoBody(state_hist1(:,ii - 1), dt, J, [0 0 0]');
    state_f2 = PropagateTwoBody(state_hist2(:,ii - 1), dt, J, [0 0 0]');
    state_f3 = PropagateTwoBody(state_hist3(:,ii - 1), dt, J, [0 0 0]');
    state_hist1(:,ii) = state_f1;
    state_hist2(:,ii) = state_f2;
    state_hist3(:,ii) = state_f3;
end

%% Plotting
close all

figure
for ii = 1:3
    subplot(3,3,3*ii - 2)
    plot(t,state_hist1(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
    title("Case 1")

    subplot(3,3,3*ii - 1)
    plot(t,state_hist2(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
    title("Case 2")

    subplot(3,3,3*ii - 0)
    plot(t,state_hist3(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
    title("Case 3")
end