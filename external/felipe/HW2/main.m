%% Problem 1

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

Jb_cg  = diag([60;90;60]);
n      = 0.5;
wb_bi0 = [0;n;0];
taub_0 = [0.001;0;0];

[tnl,wb_binl] = ode45(@(t,x) euler(t,x,Jb_cg,taub_0),[0 60],wb_bi0,options);
[tl,wb_bil]   = ode45(@(t,x) eulerlin(t,x,Jb_cg,taub_0,n),[0 60],[0;0;0],options);

figure()
tiledlayout(3,1)
nexttile
plot(tnl,wb_binl(:,1),'LineWidth',2); hold on;
plot(tl,wb_bil(:,1) + wb_bi0(1),'--','LineWidth',2)
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('$\omega^{b}_{b/i_{x}}$ (rad/s)', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);
nexttile
plot(tnl,wb_binl(:,2),'LineWidth',2); hold on;
plot(tl,wb_bil(:,2) + wb_bi0(2),'--','LineWidth',2)
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('$\omega^{b}_{b/i_{y}}$ (rad/s)', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);
nexttile
plot(tnl,wb_binl(:,3),'LineWidth',2); hold on;
plot(tl,wb_bil(:,3) + wb_bi0(3),'--','LineWidth',2)
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('$\omega^{b}_{b/i_{z}}$ (rad/s)', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);

%% Problem 2

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

G      = 6.6743e-11;                    % m3/kgs2
m_e    = 5.97219e+24;                   % kg 
Jb_cg  = diag([90 70 60]);              % kgm2
a_0    = 400e+3 + 6.378137e+6;          % m
Torbs  = 2*pi*sqrt((a_0^3)/(G*m_e));
n      = -2*pi/Torbs;
h0     = 2*(5e-5)/(n*0.5*pi/180);
gamma  = 2*Jb_cg(3,3)*n/h0;
kp     = 1.5;
tau    = 70;
kd     = kp*tau;

[t,x]  = ode45(@(t,x) mwlin(t,x,Jb_cg,n,h0,gamma,kp,kd),[0 1*Torbs],[0;0;0;0.0001;0;0.0001],options);

figure()
tiledlayout(3,1)
nexttile
plot(t,x(:,1)*180/pi,'LineWidth',2); hold on;
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('Roll (deg)', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);
nexttile
plot(t,x(:,2)*180/pi,'LineWidth',2); hold on;
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('Pitch (deg)', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);
nexttile
plot(t,x(:,3)*180/pi,'LineWidth',2); hold on;
grid on
grid minor
xlabel('time (s)', 'Interpreter','latex')
ylabel('Yaw (deg)', 'Interpreter','latex')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex','FontSize',20);

%% Functions

% Problem 1
function xdot = euler(t,x,Jb,tau)
    % Unpacking states
    wb_bi = x;
    % Propagation of dynamics
    xdot  = inv(Jb)*(-cross(wb_bi,Jb*wb_bi) + tau);
end

function xdot = eulerlin(t,x,J,tau,n)
    % Unpacking states
    w = x;
    xdot    = zeros(3,1);
    % Propagation of dynamics
    xdot(1) = (tau(1) + n*(J(2,2)-J(3,3))*w(3))/(J(1,1));
    xdot(2) = (tau(2))/(J(2,2));
    xdot(3) = (tau(3) + n*(J(1,1)-J(2,2))*w(1))/(J(3,3));
end

% Problem 2
function xdot = mwlin(t,x,J,n,h0,gamma,kp,kd)
    % Unpacking states
    o    = x(1:3);
    w    = x(4:6);
    xdot = zeros(6,1);
    tau  = zeros(3,1);

    tau(1)    = -kp*o(1) - kd*w(1) + n*h0*o(1);
    tau(3)    = gamma*kp*o(1) + gamma*kd*w(1) + h0*w(1);
    
    xdot(1:3) = w;
    xdot(4)   = ((4*n*n*(J(3,3)-J(2,2)) - n*h0)*o(1) + (n*(J(2,2)-J(3,3)-J(1,1))+h0)*w(3) + tau(1))/J(1,1);
    xdot(5)   = 3*n*n*(J(3,3)-J(1,1))*o(2)/J(2,2);
    xdot(6)   = ((n*n*(J(1,1)-J(2,2)) + n*h0)*o(3) + (n*(J(1,1)-J(2,2)+J(3,3))- h0)*w(1) + tau(3))/J(3,3);
end





