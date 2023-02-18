%% Setup
clear
close all
clc

addpath("..\util\")

%% Problem Initalization

% Rotation formulations
J1 = 60;
J2 = 90;
J3 = 60;
J = blkdiag(J1, J2, J3); % Inertia matrix
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
out = sim("simulink/h2p1_sim");

% Extract information
r_inertial_hist = squeeze(out.pos);
q_inertial2body_hist = squeeze(out.quat);
w_body_hist = squeeze(out.w);
tout = out.tout;

% Find the predicted deviations from nominal
k1 = (J2-J3)/J1;
k3 = (J1-J2)/J3;
Omega = sqrt(-1*k1*k3*w_b_0(2)^2);
domega = zeros(2,length(tout));
for ii = 1:length(tout)
    A = [-sin(Omega*tout(ii)), sqrt(-k1/k3)*cos(Omega*tout(ii)); sqrt(-k3/k1)*cos(Omega*tout(ii)), -sin(Omega*tout(ii))];
    domega(:,ii) = (1/Omega)*([0 1; 1 0] - A)*[J1^-1, 0; 0, J3^-1]*[T(1); T(3)];
end

%% Plotting

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(tout, w_body_hist(ii,:))
    if(ii == 1)
        hold on
        plot(tout,domega(1,:))
    end

    if(ii == 3)
        hold on
        plot(tout,domega(2,:))
    end
    xlabel('Time [s]')
    ylabel('Body Rate [rad/sec]')
end
