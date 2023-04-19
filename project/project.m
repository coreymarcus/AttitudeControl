%% Setup
clear
close all
clc

addpath("..\util\")

modelname = "simulink\project_sim.slx";

%% Symbolic math for some functions
perform_symbolic = true;
if(perform_symbolic)
    w_hat = sym('w_hat',[3 1],'real');
    syms dt tau 'real'

    wx = CrossProductMat(w_hat);
    A = ExpmSkewSym(-wx*dt);
    B = int(ExpmSkewSym(-wx*(dt - tau)),tau,0,dt);

    matlabFunction(A,'File','PaPaFunc.m');
    matlabFunction(B,'File','PaPdomegaFunc.m');
end

%% Problem Initalization

Sigma_a = 1E-2*eye(3); % Measurement noise for q_inertial2body
underweight = 10;

% Process noise for filter
Q_filter = blkdiag(1E-4*eye(3),1E-9*eye(3));

% Initial uncertainty for filter
Phat0 = 10*Q_filter;

% Flight software frequency
FSW_freq = 1/10;

% Gravity
mu = 398600; %km3/s2

% Orbital elements
a = 6371 + 400;
e = 0;
i = 0;
Ohm = 0;
w = 0;
theta = 0;
[r_inertial_0, v_inertial_0] = OE2State(a, e, i, Ohm, w, theta);

% Find orbit rate
n = sqrt(mu/a^3);

% Find LVLH rotation rate
w_LVLH_wrt_inertial_in_LVLH = [0, -n, 0]';

% Initial rotation between LVLH and inertial
x_LVLH_inertial_0 = v_inertial_0/norm(v_inertial_0);
z_LVLH_inertial_0 = -r_inertial_0/norm(r_inertial_0);
y_LVLH_inertial_0 = cross(z_LVLH_inertial_0,x_LVLH_inertial_0);
T_inertial2LVLH_0 = [x_LVLH_inertial_0'; y_LVLH_inertial_0'; z_LVLH_inertial_0'];
q_inertial2LVLH_0 = DCM2Quat(T_inertial2LVLH_0);

% Rotation formulations
J = 1E-4*[24181836 3783405 3898808
    3783405 37621803 -1171849
    3898808 -1171849 51576634];
w_b_LVLH_0 = [0 0 0]'; % Initial LVLH rotation rate, rad/sec
q_LVLH2body_0 = [0.028, -0.0788, 0.1141, 0.9899]'; % Initial attitude quaternion
%q_LVLH2body_f = q_LVLH2body_0;
q_LVLH2body_f = [-0.0607, -0.0343, -0.7045, 0.7062]'; % Attitude quaternion at end of the manuever

% Initial pose and rate in inertial
q_inertial2body_0 = QuatProduct(q_LVLH2body_0,q_inertial2LVLH_0);
w_body_wrt_inertial_0 = QuatTransform(q_LVLH2body_0,w_LVLH_wrt_inertial_in_LVLH) + w_b_LVLH_0;

% Final simulation time
Tf = 2*7110;
% Tf = 100;
% Tf = 10000;

% Final manuever time
Tf_man = 7110;

% CMG momentum
h0 = 4881;

% Maximum CMG rates
rate_max = Inf*(pi/180); % Rad/sec

%% Sample discrete time noise
tsample = 0:FSW_freq:Tf_man;
Nsample = length(tsample);
q_meas_noise = mvnrnd(zeros(3,1),Sigma_a./underweight,Nsample);
simin.q_meas_noise = timeseries(q_meas_noise,tsample);

%% Design the maneuver in the LVLH frame

% Change in quaternion
dq_LVLH = QuatProduct(q_LVLH2body_f,QuatInv(q_LVLH2body_0));

% Euler axis and angle change
[dtheta_LVLH, dn_LVLH] = Quat2AxisAngle(dq_LVLH);

% Find angular rate in rad/sec
w_b_LVLH_man = dtheta_LVLH/Tf_man*dn_LVLH;

%% Nonlinear controller design

kp_nonlin = 0.1;
kd_nonlin = 1;

%% Main

use_CMG = true;
out_data = sim(modelname);

%% Extract Information

w_body_inertial = out_data.w;
w_body_inertial_est = out_data.w_body_est;
q_inertial2LVLH = out_data.q_inertial2LVLH;
q_inertial2body = out_data.quat;
q_inertial2body_est = out_data.q_inertial2body_est;
CMG_rates = out_data.CMG_rates;
% err_quat = squeeze(out_data.error_quat)';
w_body_inertial_ref = out_data.ref_rate';
CMG_h = squeeze(out_data.CMG_h)';

% Squeeze for some reason
q_inertial2body_est.Data = squeeze(q_inertial2body_est.Data)';
w_body_inertial_est.Data = squeeze(w_body_inertial_est.Data)';

% % Find quaternion from body to LVLH
% q_LVLH2body = zeros(size(q_inertial2LVLH));
% for ii = 1:length(tout)
%     q_LVLH2body(:,ii) = QuatProduct(q_inertial2body(:,ii),QuatInv(q_inertial2LVLH(:,ii)));
% end

%% Plotting Part 3

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(out_data.error_quat.Time,out_data.error_quat.Data(:,ii),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel("$\delta q$","Interpreter","latex")
    grid on
    title("Estimated Pointing Error")
end

figure
for ii = 1:4
    subplot(4,1,ii)
    hold on
    plot(q_inertial2body.Time, q_inertial2body.Data(:,ii),"LineWidth",2)
    plot(q_inertial2body_est.Time, q_inertial2body_est.Data(:,ii),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel('Quat Element',"Interpreter","latex")
    legend("Actual","Estimate")
    grid on
end

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(w_body_inertial.Time, w_body_inertial.Data(:,ii),"LineWidth",2)
    plot(w_body_inertial_est.Time, w_body_inertial_est.Data(:,ii),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel('Angular Rate',"Interpreter","latex")
    legend("Actual","Estimate")
    grid on
end

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(w_body_inertial.Time, w_body_inertial_ref.Data(:,ii) - w_body_inertial.Data(:,ii),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel('$\delta \omega$ [rad/sec]',"Interpreter","latex")
    legend("Actual","Estimate")
    grid on
end

figure
for ii = 1:4
    subplot(4,2,2*ii-1)
    plot(CMG_rates.Time,CMG_rates.Data(:,ii))
    xlabel('Time [sec]',"Interpreter","latex")
    ylabel(strcat("$\dot{\alpha}$ ",num2str(ii)," [rad/sec]"),"Interpreter",'latex')
    grid on

    subplot(4,2,2*ii)
    plot(CMG_rates.Time,CMG_rates.Data(:,ii+4))
    xlabel('Time [sec]',"Interpreter","latex")
    ylabel(strcat("$\dot{\beta}$ ",num2str(ii)," [rad/sec]"),"Interpreter",'latex')
    grid on

end

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(CMG_h.Time, abs(CMG_h.Data(:,ii)),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel('$|h_{CMG}|$ [kg-m\textsuperscript{2}/sec]',"Interpreter","latex")
    grid on
end
