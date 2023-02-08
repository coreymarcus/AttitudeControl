clear
close all
clc

rng(3)

% Intialize two quaternions
q1 = rand(4,1);
q1 = q1/norm(q1);
q2 = rand(4,1);
q2 = q2/norm(q2);

% Compose quaternions
q3 = QuatProduct(q2,q1);
T3 = Quat2DCM(q3)

% Compose DCMs
T1 = Quat2DCM(q1);
T2 = Quat2DCM(q2);
T3_2 = T2*T1