%% Setup
clear
close all
clc

addpath("..\util\")

%% Problem Initalization

createmovie = true;

% Rotation formulations
J = blkdiag(100, 75, 50); % Inertia matrix
w_b_0_1 = [0 .1 0]'; % Initial rotation rate, rad/sec
w_b_0_2 = [0.01 .0001 0.0001]'; % Initial rotation rate, rad/sec
w_b_0_3 = [0.0001 .01 0.0001]'; % Initial rotation rate, rad/sec
q_inertial2body_0 = [0 0 0 1]'; % Initial attitude quaternion

% We don't actually need position and velocity but easier to interface with
% the code this way
a = 6371 + 500;
e = 0;
i = 0;
Ohm = 0;
w = 0;
theta = 0;

% Final time
Tf = 2.5*60*60;
% Tf = 25;

%% Main

% Initialize x y z unit vectors
xx = [1 0 0]';
yy = [0 1 0]';
zz = [0 0 1]';

% Convert OE
[r_inertial, v_inertial] = OE2State(a, e, i, Ohm, w, theta);

% Propagation timestep
dt = 10;
% dt = 1;

% Time vector
t = 0:dt:Tf;
Nt = length(t);

% Initialize state history vector
state_hist1 = zeros(13,Nt);
state_hist1(1:3,1) = r_inertial;
state_hist1(4:6,1) = v_inertial;
state_hist1(7:10,1) = q_inertial2body_0;
state_hist2 = state_hist1;
state_hist3 = state_hist1;
state_hist1(11:13,1) = w_b_0_1;
state_hist2(11:13,1) = w_b_0_2;
state_hist3(11:13,1) = w_b_0_3;

% Initialize triad history
xx_inertial1 = zeros(3,Nt);
xx_inertial2 = zeros(3,Nt);
xx_inertial3 = zeros(3,Nt);
yy_inertial1 = zeros(3,Nt);
yy_inertial2 = zeros(3,Nt);
yy_inertial3 = zeros(3,Nt);
zz_inertial1 = zeros(3,Nt);
zz_inertial2 = zeros(3,Nt);
zz_inertial3 = zeros(3,Nt);

% Initialize movie frames
M1(Nt-1) = struct('cdata',[],'colormap',[]);
M2(Nt-1) = struct('cdata',[],'colormap',[]);
M3(Nt-1) = struct('cdata',[],'colormap',[]);

figure;

% Propagate
for ii = 2:Nt
    state_f1 = PropagateTwoBody(state_hist1(:,ii - 1), dt, J, [0 0 0]');
    state_f2 = PropagateTwoBody(state_hist2(:,ii - 1), dt, J, [0 0 0]');
    state_f3 = PropagateTwoBody(state_hist3(:,ii - 1), dt, J, [0 0 0]');
    state_hist1(:,ii) = state_f1;
    state_hist2(:,ii) = state_f2;
    state_hist3(:,ii) = state_f3;

    % Rotate triads
    xx_inertial1(:,ii) = QuatTransform(QuatInv(state_hist1(7:10,ii)),xx);
    xx_inertial2(:,ii) = QuatTransform(QuatInv(state_hist2(7:10,ii)),xx);
    xx_inertial3(:,ii) = QuatTransform(QuatInv(state_hist3(7:10,ii)),xx);
    yy_inertial1(:,ii) = QuatTransform(QuatInv(state_hist1(7:10,ii)),yy);
    yy_inertial2(:,ii) = QuatTransform(QuatInv(state_hist2(7:10,ii)),yy);
    yy_inertial3(:,ii) = QuatTransform(QuatInv(state_hist3(7:10,ii)),yy);
    zz_inertial1(:,ii) = QuatTransform(QuatInv(state_hist1(7:10,ii)),zz);
    zz_inertial2(:,ii) = QuatTransform(QuatInv(state_hist2(7:10,ii)),zz);
    zz_inertial3(:,ii) = QuatTransform(QuatInv(state_hist3(7:10,ii)),zz);

    % Create movie frame
    if(createmovie)
        clf
        hold on
        quiver3(0,0,0,xx_inertial1(1,ii),xx_inertial1(2,ii),xx_inertial1(3,ii))
        quiver3(0,0,0,yy_inertial1(1,ii),yy_inertial1(2,ii),yy_inertial1(3,ii))
        quiver3(0,0,0,zz_inertial1(1,ii),zz_inertial1(2,ii),zz_inertial1(3,ii))
        grid on
        legend("Body X","Body Y","Body Z")
        axis([-1 1 -1 1 -1 1])
        title("Case 1")
        view(45,45)
        M1(ii-1) = getframe;

        clf
        hold on
        quiver3(0,0,0,xx_inertial2(1,ii),xx_inertial2(2,ii),xx_inertial2(3,ii))
        quiver3(0,0,0,yy_inertial2(1,ii),yy_inertial2(2,ii),yy_inertial2(3,ii))
        quiver3(0,0,0,zz_inertial2(1,ii),zz_inertial2(2,ii),zz_inertial2(3,ii))
        grid on
        legend("Body X","Body Y","Body Z")
        axis([-1 1 -1 1 -1 1])
        title("Case 2")
        view(45,45)
        M2(ii-1) = getframe;

        clf
        hold on
        quiver3(0,0,0,xx_inertial3(1,ii),xx_inertial3(2,ii),xx_inertial3(3,ii))
        quiver3(0,0,0,yy_inertial3(1,ii),yy_inertial3(2,ii),yy_inertial3(3,ii))
        quiver3(0,0,0,zz_inertial3(1,ii),zz_inertial3(2,ii),zz_inertial3(3,ii))
        grid on
        legend("Body X","Body Y","Body Z")
        axis([-1 1 -1 1 -1 1])
        view(45,45)
        title("Case 3")
        M3(ii-1) = getframe;
    end
end

%% Movies
if(createmovie)
    %moviefig = figure;
    %movie(moviefig, M1)
    %movie(moviefig, M2)
    %movie(moviefig, M3)

    v1 = VideoWriter("moive1","Motion JPEG AVI");
    v1.FrameRate = 10;
    open(v1);
    writeVideo(v1,M1);
    close(v1);
    v2 = VideoWriter("movie2","Motion JPEG AVI");
    v2.FrameRate = 10;
    open(v2);
    writeVideo(v2,M2);
    close(v2);
    v3 = VideoWriter("movie3","Motion JPEG AVI");
    v3.FrameRate = 10;
    open(v3);
    writeVideo(v3,M3);
    close(v3);
end

%% Plotting
close all

figure
for ii = 1:3
    subplot(3,3,3*ii - 2)
    plot(t,state_hist1(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
    title("Case 1")

    subplot(3,3,3*ii - 1)
    plot(t,state_hist2(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
    title("Case 2")

    subplot(3,3,3*ii - 0)
    plot(t,state_hist3(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
    title("Case 3")
end