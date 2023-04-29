% This script finds derivatives symbolically

clear
close all
clc

addpath("..\util\")

%% Declare variables

q0 = sym("q",[4 1],"real"); % Nominal inertial to body quaternion
a = sym("a",[3 1],"real"); % Nominal error peturbation
r = sym("r",[3 1],"real"); % True inertial position

%% Model measurement

% Magnetic field in inertial
B_inertial = DipoleMagneticField(r);

% Convert a to a quaternion
dq = [0.5*a; 1];
dq = dq/norm(dq);

% B in body
B_body = NativeQuatTransform(QuatProduct(dq,q0),B_inertial);

% Derivative
pBpa = jacobian(B_body,a);

% Convert to functions
matlabFunction(pBpa,"File","partialBwrtA","Optimize",true);