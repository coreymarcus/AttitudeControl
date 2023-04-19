clear
close all
clc

addpath("..\util\")

% Parameters
dt = 0.1;
q0 = rand(4,1);
q0 = q0/norm(q0);
w0 = rand(3,1);
J = rand(3,3) + 3*eye(3);
J = J + J'


% Propagate with ODE45
tic
[qf, wf] = AttitudePropagate(q0, w0, J, dt, "ode45")
toc

% Propagate with RK4
tic
[qf, wf] = AttitudePropagate(q0, w0, J, dt, "RK4")
toc
