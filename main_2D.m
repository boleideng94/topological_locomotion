%% clear all
clear all
global Nx Ny a ks kt kb kl lrb mass J
global damp_theta damp_u Time thst mass_act
global Time_pre Time_wait Time_switch cycle
global actuator_tube progress lp la
global friction_dt friction_time friction_tube label
global cycle_number cycle_number_current wheel_angle
global mu_0 mu_1 pin_support mu_pin

%% experiment type
pin_support = 1; % choose 0 if it is locomotion test (wheels supported)
                   % choose 1 if it is wave propagation test (pin supported)
                   
%% unti-cell parameters
% hinges
lh_inter = 6.*1e-3; % inter-hinge length
lh_intra = 3.5e-3; % intra-hinge length
ws = 1.6e-2; % width of the shims
d_inter = 0.127e-3; % thikness of the inter-hinge
d_intra = 0.25e-3; % thickness of the intra-hinge
E_shim = 4.3e9; % stiffness of the shim
I_inter = ws*d_inter^3/12;
I_intra = ws*d_intra^3/12;

% unit dimensions
la = 19e-3; % arm length of the cross
lp = 32.5e-3+lh_intra; % post distance

% rubber bands
lrb = 27.5*1e-3; % rubber band initial length
Erb = 7e6; % stiffness of rubber band
wrb = 2*0.3e-3; % width of the rubber band
d_rb = 2e-3; % thickness of the rubber band

% spring stiffness calculated from geometrical parameters
ks = 12*E_shim*I_inter/lh_inter^3; % shearing spring
kt = 0.0103; % torsional spring
kb = Erb*wrb*d_rb/lrb; % linear spring for rubber band 
kl = 100e3; % longitudinal spring

% other basic parameters
mass = 31.5e-3; % mass of the unit cell 
mass_act = 35e-3; % mass of the actuator
J = 1/3*mass*(19e-3)^2*10; % rotational inertia
damp_theta = 0; % viscous damping
damp_u = 0.0; % viscous damping

% calculated values
thst = (2-2*lrb/lp-8*kt/(lp^2*kb))^0.5;
disp(['The static angle of the unit cell, theta_st = ',num2str(thst*180/pi),' deg']);

%% structure and actuation parameters

Nx = 5; % number of unit cells in x-direction N_x
Ny = 1;  % number of unit cells in y-direction N_y
actuator_tube = [1,1]; % position of the actuator(s)
wheel_angle = (0).*ones([Ny,Nx]); % wheel principal angle (default as 0)

% simulation time & animation play rate
Time_pre = 0.2; % prepare time (to settle down transient vibrations)
Time_switch = 0.5*Nx; % actuation time of the actuator
Time_wait = 0.5*Nx; % waiting time at the end of actuation
cycle = 2; % number of half cycles
Time = Time_pre + (Time_switch + Time_wait)*cycle*length(actuator_tube(:,1)); % total time

% flags and trackers
cycle_number = 0; progress = 0; cycle_number_current = 0;
label = 1; flag_cycle = 1;

% friction parameters
friction_dt = 1e-5; friction_time = 0; friction_tube = zeros([Ny,Nx,3]);

% friction of the wheel
% \mu_s = 0.95 + 7.8v_s (Eq. S40)
mu_0 = 0.95; mu_1 = 7.8;

% friction of the pin 
mu_pin = 0.15;

%% initial condition and ODE solver
y0 = 0e-5.*rand(8*Nx*Ny,1);
for i = 1:Ny
    for j = 1:Nx
        posi = ((i-1)*Nx + j-1)*8;
        Th_cen=posi+3;
        y0(Th_cen) = thst; % set initial angle to theat_st
    end
end

tic;
[T,Yout] = ode45('ODE_2D',[0,Time], y0); % run ode45 solver
[Timelength, Distance]=size(Yout); % output data
timeelapse=toc;
disp(['Elapsed time = ',num2str(timeelapse),' s']);

%% Data processing

FrameN = 100; % totoal frame 
pre_len =  round(interp1(T,1:Timelength,Time_pre));
clear DispX DispY thetaR theta_R UX UY beta_R betaR
for i = 1:Ny
    for j = 1:Nx
        posi = ((i-1)*Nx + j-1)*8;
        DispX(i,j,:) = Yout(:,posi+1);                                      % mark here
        DispY(i,j,:) = Yout(:,posi+2);
        theta_R(i,j,:) = Yout(:,posi+3);
        beta_R(i,j,:) = Yout(:,posi+4);
    end
end
count = 0;

% find array of time point with total number of FrameN
timepoints = round(interp1(T(pre_len:end),pre_len:length(T),linspace(T(pre_len),T(end-1),FrameN)));

for t=1:FrameN
    count = count + 1;
    for i=1:Ny
        for j=1:Nx
            UX(i,j,t) = DispX(i,j,timepoints(t));%x diplacement for the center
            UY(i,j,t) = DispY(i,j,timepoints(t));%y displacment for the center
            thetaR(i,j,t) = theta_R(i,j,timepoints(t));%internal rotation angle
            betaR(i,j,t) = beta_R(i,j,timepoints(t));%global rotation angle
        end
    end
end
Tout = T(timepoints);
step_len = (UX(1,1,1) - UX(1,1,end))/cycle; % step length

% save data
if pin_support == 1
    folder = 'wave_propagation';
else
    folder = 'locomotion';
end
mkdir(['data/',folder]);
file = ['Nx',num2str(Nx),'_Ny',num2str(Ny),'_inter_hinge',num2str(lh_inter*1e3)];
save(['data/',folder,'/',file,'.mat'],'UX','UY','thetaR','betaR','lrb','lh_inter','lh_intra','ks','kb','kt','kl','E_shim','ws',...
    'd_inter','d_intra','a','Erb','wrb','d_rb','Nx','Ny','FrameN','Tout','actuator_tube','step_len','mass','mass_act');
close all

%% movies

mkdir(['movies/',folder]);
fig = figure(1);
set(fig,'position', [0, 0,(Nx+0.5)*200,(Ny+0.5)*200]);
set(gcf,'Color','w');
myVideo = VideoWriter(['movies/',folder,'/',file],'MPEG-4');
myVideo.FrameRate = 30;
open(myVideo);

a = la*4+0.004;
for i=1:Ny
    for j=1:Nx
        X0(i,j) = j*a*cos(thst)*1.1;
        Y0(i,j) = i*a*cos(thst)*1.1;
    end
end

for t = 1:FrameN 
    clf
    for i = 1:Ny
        for j = 1:Nx
            center = [X0(i,j),Y0(i,j)] + [UX(i,j,t),UY(i,j,t)];
            theta = thetaR(i,j,t);
            beta = betaR(i,j,t);
            element_plotter(theta,beta,center,a);
            if sum((i == actuator_tube(:,1)).*(j == actuator_tube(:,2)))
                tt = linspace(0,2*pi,100);
                x = center(1)+a/14.*cos(tt); y = center(2)+a/14.*sin(tt); 
                fill(x,y,'M','LineWidth',1);
            end
        end
    end
    axis equal
    axis off
    xlim([-2*a,(Nx+1)*a]); ylim([-0.5*a,(Ny+1)*a]);
    set(gca,'LooseInset',[0,0.,0,0]);
    frame = getframe(gcf);
    writeVideo(myVideo,frame);
    pause(0.00001);
end
close(myVideo);

