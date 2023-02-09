%% Setup

clear
close all
clc

addpath("..\util\")

%% Problem Initalization

createmovie = true;

% Rotation formulations
J = blkdiag(100, 100, 60); % Inertia matrix
w_b_0 = pi*[1 0 10]'/180; % Initial rotation rate, rad/sec

% Inertial z axis in body frame
i_z_in_b = J*w_b_0;
i_z_in_b = i_z_in_b/norm(i_z_in_b);
i_y_in_b = [0 1 0]';

% We don't actually need position and velocity but easier to interface with
% the code this way
a = 6371 + 500;
e = 0;
i = 0;
Ohm = 0;
w = 0;
theta = 0;

% Final time
Tf = 3000;

%% Main

% Find missing axis
i_x_in_b = cross(i_y_in_b,i_z_in_b);

% Setup initial transformation matrix
T_inertial2body0 = [i_x_in_b';
    i_y_in_b';
    i_z_in_b'];
q_inertial2body_0 = DCM2Quat(T_inertial2body0);

% Initialize z unit vectors
zz = [0 0 1]';

% Convert OE
[r_inertial, v_inertial] = OE2State(a, e, i, Ohm, w, theta);

% Propagation timestep
% dt = 10;
dt = 1;

% Time vector
t = 0:dt:Tf;
Nt = length(t);

% Initialize state history vector
state_hist = zeros(13,Nt);
state_hist(1:3,1) = r_inertial;
state_hist(4:6,1) = v_inertial;
state_hist(7:10,1) = q_inertial2body_0;
state_hist(11:13,1) = w_b_0;

% Initialize zz history
zz_inertial = zeros(3,Nt-1);

% Initialize movie frames
M(Nt-1) = struct('cdata',[],'colormap',[]);
figure;

% Propagate
for ii = 2:Nt
    state_f = PropagateTwoBody(state_hist(:,ii - 1), dt, J, [0 0 0]');
    state_hist(:,ii) = state_f;

    % Rotate vector
    zz_inertial(:,ii-1) = QuatTransform(QuatInv(state_hist(7:10,ii)),zz);

    % Create movie frame
    if(createmovie)
        clf
        hold on
        quiver3(0,0,0,zz_inertial(1,ii-1),zz_inertial(2,ii-1),zz_inertial(3,ii-1))
        grid on
        legend("Body Z")
        axis([-1 1 -1 1 -1 1])
        view(45,45)
        M(ii-1) = getframe;
    end
end

%% Plotting
figure
for ii = 1:3
    subplot(3,1,ii)
    plot(t,state_hist(10+ii,:))
    xlabel("Time [sec]")
    ylabel("Body Rate [rad/sec]")
end

if(createmovie)
    v1 = VideoWriter("movie_p3","Motion JPEG AVI");
    v1.FrameRate = 100;
    open(v1);
    writeVideo(v1,M);
end