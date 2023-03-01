%% Run sim

clc; clear;

G      = 6.6743e-11;                    % m3/kgs2
m_e    = 5.97219e+24;                   % kg 
a_0    = 400e+3 + 6.378137e+6;          % m
i_0    = 0;                             % rad
rvi_0  = ogoe2rv([a_0, 0, i_0, 0, 0, 3*pi/2],G*m_e);
Torbs  = 2*pi*sqrt((a_0^3)/(G*m_e));
worbs  = 2*pi/Torbs;
worbs  = worbs*180/pi;
Jb_cg  = diag([90 70 60]);              % kgm2
acc    = 0.5*pi/180;
tf     = 5*Torbs;

bi_x0  = rvi_0(4:6)/norm(rvi_0(4:6));
bi_y0  = -rvi_0(1:3)/norm(rvi_0(1:3));
bi_z0  = cross(bi_x0,bi_y0);
Tb_i0  = [bi_x0 bi_y0 bi_z0]';
qb_i0  = rotm2quat(Tb_i0);
wb_bi0 = [0.0001 0.0001 worbs]'*pi/180; % rad/s
x0     = [rvi_0;qb_i0;wb_bi0]';
angle0 = [0 0 0]';

out    = sim('spacecraft2');

%% Figures
qerror = out.error(:,1:4);
t      = out.tout;

dtheta1 = 2*atan2(qerror(:,1),qerror(:,4))*180/pi;
dtheta2 = 2*atan2(qerror(:,2),qerror(:,4))*180/pi;
dtheta3 = 2*atan2(qerror(:,3),qerror(:,4))*180/pi;

tiledlayout(3,1)
nexttile
plot(t,dtheta1,'LineWidth',2)
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('$\delta \theta_1 (^\circ)$', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);
nexttile
plot(t,dtheta2,'LineWidth',2)
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('$\delta \theta_2 (^\circ)$', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);
nexttile
plot(t,dtheta3,'LineWidth',2)
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('$\delta \theta_3 (^\circ)$', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);
