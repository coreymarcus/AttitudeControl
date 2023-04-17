%% Setup
clear
close all
clc

addpath("..\util\")

getfullname(Simulink.findBlocks('midterm_sim_R2021a', 'MaskType', 'Replaced Block'))

%% Problem Initalization

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
J = [24181836 3783405 3898808
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

%% Design the maneuver in the LVLH frame

% Change in quaternion
dq_LVLH = QuatProduct(q_LVLH2body_f,QuatInv(q_LVLH2body_0));

% Euler axis and angle change
[dtheta_LVLH, dn_LVLH] = Quat2AxisAngle(dq_LVLH);

% Find angular rate in rad/sec
w_b_LVLH_man = dtheta_LVLH/Tf_man*dn_LVLH;

%% Nonlinear controller design

kp_nonlin = 500;
kd_nonlin = 500000;

%% Main

use_CMG = false;
out_no_CMG = sim("simulink\midterm_sim.slx");

use_CMG = true;
out_w_CMG = sim("simulink\midterm_sim.slx");

%% Extract information for Part 1

tout_no_CMG = out_no_CMG.tout;
q_LVLH2body_ref = squeeze(out_no_CMG.ref_quat_LVLH);
w_body_LVLH_ref = squeeze(out_no_CMG.ref_rate_LVLH);
q_inertial2body_ref = squeeze(out_no_CMG.ref_quat);
w_body_inertial_ref_no_CMG = squeeze(out_no_CMG.ref_rate);


%% Extract Informatin for Part 2

err_quat_no_CMG = squeeze(out_no_CMG.error_quat);
q_inertial2body_no_CMG = squeeze(out_no_CMG.quat);
w_body_inertial_no_CMG = squeeze(out_no_CMG.w);
q_inertial2LVLH_no_CMG = squeeze(out_no_CMG.q_inertial2LVLH);

% Find quaternion from body to LVLH
q_LVLH2body_no_CMG = zeros(size(q_inertial2LVLH_no_CMG));
for ii = 1:length(tout_no_CMG)
    q_LVLH2body_no_CMG(:,ii) = QuatProduct(q_inertial2body_no_CMG(:,ii),QuatInv(q_inertial2LVLH_no_CMG(:,ii)));
end

%% Extract Information for Part 3

tout_w_CMG = out_w_CMG.tout;
w_body_inertial_w_CMG = squeeze(out_w_CMG.w);
q_inertial2LVLH_w_CMG = squeeze(out_w_CMG.q_inertial2LVLH);
q_inertial2body_w_CMG = squeeze(out_w_CMG.quat);
CMG_rates = squeeze(out_w_CMG.CMG_rates);
err_quat_w_CMG = squeeze(out_w_CMG.error_quat);
w_body_inertial_ref_w_CMG = squeeze(out_w_CMG.ref_rate);
CMG_h = squeeze(out_w_CMG.CMG_h);

% Find quaternion from body to LVLH
q_LVLH2body_w_CMG = zeros(size(q_inertial2LVLH_no_CMG));
for ii = 1:length(tout_w_CMG)
    q_LVLH2body_w_CMG(:,ii) = QuatProduct(q_inertial2body_w_CMG(:,ii),QuatInv(q_inertial2LVLH_w_CMG(:,ii)));
end


%% Plotting Part 1

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(tout_no_CMG, w_body_LVLH_ref(ii,:),'LineWidth',2)
    xlabel('Time [s]','Interpreter','latex')
    ylabel('$\omega_{LVLH}$ [rad/sec]','Interpreter','latex')
    grid on
end
saveas(gcf,"latex/figs/P1Q1.pdf")

figure
for ii = 1:4
    subplot(4,1,ii)
    hold on
    plot(tout_no_CMG, q_LVLH2body_ref(ii,:),'LineWidth',2)
    xlabel('Time [s]','Interpreter','latex')
    ylabel("$q_{b \leftarrow LVLH}$","Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P1Q2.pdf")

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(tout_no_CMG, w_body_inertial_ref_no_CMG(ii,:),'LineWidth',2)
    xlabel('Time [s]','Interpreter','latex')
    ylabel('$\omega_{inertial}$ [rad/sec]','Interpreter','latex')
    grid on
end
saveas(gcf,"latex/figs/P1Q3.pdf")

figure
for ii = 1:4
    subplot(4,1,ii)
    hold on
    plot(tout_no_CMG, q_inertial2body_ref(ii,:),'LineWidth',2)
    xlabel('Time [s]')
    ylabel("$q_{b \leftarrow inertial}$",'Interpreter','latex')
    grid on
end
saveas(gcf,"latex/figs/P1Q4.pdf")

%% Plotting Part 2

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(tout_no_CMG, err_quat_no_CMG(ii,:),'LineWidth',2)
    xlabel('Time [s]','Interpreter','latex')
    ylabel("$q_{b \leftarrow \bar{b}}$","Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P2Q1.pdf")

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(tout_no_CMG, w_body_inertial_ref_no_CMG(ii,:) - w_body_inertial_no_CMG(ii,:),'LineWidth',2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel('$\delta \omega$ [rad/sec]',"Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P2Q2.pdf")

figure
for ii = 1:4
    subplot(4,1,ii)
    hold on
    plot(tout_no_CMG, q_inertial2body_no_CMG(ii,:),'LineWidth',2)
    xlabel('Time [s]','Interpreter','latex')
    ylabel("$q_{b \leftarrow inertial}$","Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P2Q3.pdf")

figure
for ii = 1:4
    subplot(4,1,ii)
    hold on
    plot(tout_no_CMG, q_LVLH2body_no_CMG(ii,:),"LineWidth",2)
    xlabel('Time [s]', "Interpreter","latex")
    ylabel("$q_{b \leftarrow LVLH}$","Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P2Q4.pdf")

%% Plotting Part 3

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(tout_w_CMG, err_quat_w_CMG(ii,:),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel("$\delta q$","Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P3Q1.pdf")

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(tout_w_CMG, w_body_inertial_ref_w_CMG(ii,:) - w_body_inertial_w_CMG(ii,:),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel('$\delta \omega$ [rad/sec]',"Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P3Q2.pdf")

figure
for ii = 1:4
    subplot(4,2,2*ii-1)
    plot(tout_w_CMG,CMG_rates(ii,:))
    xlabel('Time [sec]',"Interpreter","latex")
    ylabel(strcat("$\dot{\alpha}$ ",num2str(ii)," [rad/sec]"),"Interpreter",'latex')
    grid on

    subplot(4,2,2*ii)
    plot(tout_w_CMG,CMG_rates(ii+4,:))
    xlabel('Time [sec]',"Interpreter","latex")
    ylabel(strcat("$\dot{\beta}$ ",num2str(ii)," [rad/sec]"),"Interpreter",'latex')
    grid on

end
saveas(gcf,"latex/figs/P3Q3.pdf")

figure
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(tout_w_CMG, abs(CMG_h(ii,:)),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel('$|h_{CMG}|$ [kg-m\textsuperscript{2}/sec]',"Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P3Q4.pdf")

figure
for ii = 1:4
    subplot(4,1,ii)
    hold on
    plot(tout_w_CMG, q_LVLH2body_w_CMG(ii,:),"LineWidth",2)
    xlabel('Time [s]',"Interpreter","latex")
    ylabel("$q_{b \leftarrow LVLH}$","Interpreter","latex")
    grid on
end
saveas(gcf,"latex/figs/P3Q5.pdf")