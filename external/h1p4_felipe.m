%% Setup - Corey
clear
close all
% clc


%% Main Body
% Problem 4
global G m_e
G      = 6.6743e-11;                    % m3/kgs2
m_e    = 5.97219e+24;                   % kg
a_0    = 400e+3 + 6.378137e+6;          % m
i_0    = 0;                             % rad
rvi_0  = ogoe2rv([a_0, 0, i_0, 0, 0, 3*pi/2],G*m_e);
Torbs  = 2*pi*sqrt((a_0^3)/(G*m_e));
worbs  = 2*pi/Torbs;
worbs  = worbs*180/pi;
Jb_cg  = diag([90 70 60]);              % kgm2

bi_x0  = rvi_0(4:6)/norm(rvi_0(4:6));
bi_y0  = -rvi_0(1:3)/norm(rvi_0(1:3));
bi_z0  = cross(bi_x0,bi_y0);
Tb_i0  = [bi_x0 bi_y0 bi_z0]';
% qb_i0  = rotm2quat(Tb_i0);
qb_i0 = [.5 .5 .5 1]';
qb_i0 = qb_i0/norm(qb_i0);
wb_bi0 = [0.0001 0.0001 worbs]'*pi/180; % rad/s
x0     = [rvi_0;qb_i0;wb_bi0]';

ts      = 0:100:5*Torbs;
xs      = zeros(length(ts),length(x0));
xs(1,:) = x0;

for i = 1:length(ts)-1
    tm    = ts(i);
    tp    = ts(i+1);
    x     = xs(i,:);

    [~,x] = ode45(@(t,x) time_propagation4(t,x,Jb_cg,true), [tm tp], x);

    ri    = x(end,1:3);
    vi    = x(end,4:6);
    qb_i  = x(end,7:10);
    qb_i  = qb_i/norm(qb_i);
    wb_bi = x(end,11:13);

    xs(i+1,:) = [ri vi qb_i wb_bi];
end

figure
for ii = 1:3
    subplot(3,1,ii)
    plot(xs(:,6+ii))
end


%% Helpers
function xdot = time_propagation4(t,x,Jb,grav)
G      = 6.6743e-11;
m_e    = 5.97219e+24;
mu = G*m_e;

% Unpacking states
ri          = x(1:3);
vi          = x(4:6);
qb_i        = x(7:10);
wb_bi       = x(11:13);
xdot        = zeros(13,1);
% Propagation of position & velocity
xdot(1:3)   = vi;
xdot(4:6)   = -G*m_e*ri/(norm(ri)^3);
% Propagation of kinematics
xdot(7:10)  = schuster_mult([wb_bi; 0],qb_i)/2;
% Propagation of dynamics
if grav == true
    rb = schuster_mult(qb_i,schuster_mult([ri;0],star(qb_i)));
    rb = -rb(1:3);
    ub = 3*mu*cross(rb,Jb*rb)/(norm(rb)^5)
else
    ub = 0;
end
xdot(11:13) = inv(Jb)*(-cross(wb_bi,Jb*wb_bi) + ub);
end