% Symbolic Manipulations for this problem

clear
close all
clc

addpath("..\util\")

%% Initialize variables

domega = sym('domega',[3 1], 'real');
dtheta = sym('dtheta',[3 1], 'real');
omegabar = sym('omegabar',[3 1],'real');

%% Evolution of delta theta dot

% Error quaternion
dq = [dtheta/2; 1];

% Evolution of dtheta
dthetadot = 2*( 0.5*QuatProduct([omegabar; 0],dq) - 0.5*QuatProduct(dq,[omegabar; 0]) + 0.5*QuatProduct([domega; 0],dq) );

% Replace second order effects with zero
dthetadot = subs(dthetadot,cross(dtheta/2,domega),[0; 0; 0]);

% Care about vector component only
dthetadot = dthetadot(1:3);

% Expected output
dthetadot_exp = -cross(omegabar,dtheta) + domega;

% Compare with expected output
% pretty(dthetadot)
% pretty(simplify(dthetadot - dthetadot_exp))

%% Gravity gradient torque

% Mean motion and angular momentum
syms n J1 J2 J3 'real'
J = blkdiag(J1,J2,J3);

% down in nominal body frame
db_nom = [0 1 0]';

% down in actual body frame
db = (eye(3) - CrossProductMat(dtheta))*db_nom;

% Torque
tau_gg = simplify(3*n^2 * cross(db,J*db));

% Neglect the dtheta1*dtheta3 term
tau_gg = subs(tau_gg,dtheta(1)*dtheta(3),0);

%% Momentum wheel

% wheel inertia
syms h0 'real'
h_wheel_wrt_body = [0 0 h0]';

% Nominal rotation rate
omegabar_prob = [0 0 n]';

% Torque on spacecraft provided by wheel
tau_wheel = -cross(omegabar_prob+domega,h_wheel_wrt_body);

%% Putting it all together

% Control input
tau_cont = sym('tau_cont',[3 1],'real');

% Symbolic representation for dthetadot
dthetadotsym = sym('dthetadot',[3 1],'real');

% Change in omega dot
domegadot = inv(J)*(CrossProductMat(J*omegabar_prob) - CrossProductMat(omegabar_prob)*J)*domega + inv(J)*tau_gg + inv(J)*tau_wheel + inv(J)*tau_cont;

% We know the following
dthetadotdot = [n*dthetadotsym(2); -n*dthetadotsym(1); 0] + domegadot;

% Make some substitutions
dthetadotdot = subs(dthetadotdot,domega(1),dthetadotsym(1) - n*dtheta(2));
dthetadotdot = subs(dthetadotdot,domega(2),dthetadotsym(2) + n*dtheta(1));

%% Look at pitch and yaw controllers

% Assume h0 >> n J_max
dthetadotdot_py = dthetadotdot(1:2);
dthetadotdot_py = subs(dthetadotdot_py,n*J(2,2),0);
dthetadotdot_py = subs(dthetadotdot_py,n*J(3,3),0);
dthetadotdot_py = subs(dthetadotdot_py,n^2,0); % Proxy for n^2*J(2,2) or n^2*(J(3,3))

% Add a disturbance torque
tau_dist = sym('tau_dist',[2 1],'real');
dthetadotdot_py = dthetadotdot_py + tau_dist./[J(1,1); J(2,2)];

% Assume forms for control inputs
syms kp1 kp2 kd1 kd2 'real'
% Renato's forms
% dthetadotdot_py = subs(dthetadotdot_py,tau_cont(1), -kp1*dtheta(1) - kd1*dthetadotsym(1) + n*h0*dtheta(1));
% dthetadotdot_py = subs(dthetadotdot_py,tau_cont(2), -kp2*dtheta(2) - kd2*dthetadotsym(2) + n*h0*dtheta(2));

% My hacky cancellation forms
dthetadotdot_py = subs(dthetadotdot_py,tau_cont(1), -kp1*dtheta(1) - kd1*dthetadotsym(1) + n*h0*dtheta(1) + h0*dthetadotsym(2));
dthetadotdot_py = subs(dthetadotdot_py,tau_cont(2), -kp2*dtheta(2) - kd2*dthetadotsym(2) + n*h0*dtheta(2) - h0*dthetadotsym(1));

% % My hacky cancellation forms
% dthetadotdot_py = subs(dthetadotdot_py,tau_cont(1), -kp1*dtheta(1) - kd1*dthetadotsym(1) + h0*dthetadotsym(2));
% dthetadotdot_py = subs(dthetadotdot_py,tau_cont(2), -kp2*dtheta(2) - kd2*dthetadotsym(2) - h0*dthetadotsym(1));

% Perform a laplace transform
syms s 'real'
dtheta_s = sym('dtheta_s',[2 1],'real');
dthetadotdot_py_s = subs(dthetadotdot_py,dtheta(1), dtheta_s(1));
dthetadotdot_py_s = subs(dthetadotdot_py_s,dtheta(2), dtheta_s(2));
dthetadotdot_py_s = subs(dthetadotdot_py_s,dthetadotsym(1), s*dtheta_s(1));
dthetadotdot_py_s = subs(dthetadotdot_py_s,dthetadotsym(2), s*dtheta_s(2));

% Find A such that A*dtheta(s) = tau_dist(s)
pretty(collect(simplify(s^2*dtheta_s - dthetadotdot_py_s), dtheta_s))

% Inspection of the above helps us to find the following components of A
% We also manually remove extraneous J1 and J2 terms
A11 = kp1 + kd1*s + J(1,1)*s^2
A12 = -J(1,1)*n*s
A21 = J(2,2)*n*s + J(1,1)*n*s
A22 = kp2 + kd2*s + J(2,2)*s^2 - J(1,1)*n^2

% Assume h0 >> nJ_max
A12_sub = subs(A12,J(1,1)*n*s,0);
A21_sub = subs(A21,J(2,2)*n*s,0);
A21_sub = subs(A21_sub,J(1,1)*n*s,0);

% The detirminant will tell us the poles of the system
% A = [A11, A12_sub,;
%     A21_sub, A22];
A = [A11, A12,;
    A21, A22];
detA = det(A);

% Collect terms of s
T1 = collect(simplify(detA), s);
C1 = coeffs(T1,s);

% We would like our system to be almost decoupled, in that case
detA_decouple = A11*A22;
T2 = collect(simplify(detA_decouple), s);
C2 = coeffs(T2,s);

% Our goal is to minimize the difference between C1 and C2
diffC = simplify(C1 - C2)

% We notice that the poles are decoupled if
h01 = 0.5*(2*J(1,1)*n + J(2,2)*n + sqrt((-2*J(1,1)*n - J(2,2)*n)^2 - 4*J(1,1)^2));
h02 = 0.5*(2*J(1,1)*n + J(2,2)*n + sqrt((-2*J(1,1)*n - J(2,2)*n)^2 - 4*J(1,1)^2));

% Make substitutions