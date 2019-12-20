%% parameters initialization tool
% random variables parameter:
clear
clc
close all

%% add paths for functions
% addpath('Partials');
% addpath('Model');
% addpath('Validation');
% addpath('Matching_JCBB');

rng(124)

Nmax=6000;
max_landmarks=400;

times_history=struct('iter',zeros(Nmax,5));

simulator_history=struct('state',zeros(24,Nmax));

simulator_history.measures=cell(Nmax,1);

state_history=struct('xk',zeros(24,Nmax),'xk1',zeros(24,Nmax));

visibility_check.filter_meas=cell(Nmax,1);
visibility_check.H_J=cell(Nmax,1);
visibility_check.Z_J=cell(Nmax,1);
visibility_check.vis=cell(Nmax,1);
visibility_check.match=cell(Nmax,1);

landmarks_map.Sk=cell(Nmax,1);
landmarks_map.Sk1=cell(Nmax,1);
landmarks_map.new=cell(Nmax,1);

covariance_mat.Prrk=cell(Nmax,1);
covariance_mat.Prmk=cell(Nmax,1);
covariance_mat.Pmmk=cell(Nmax,1);

covariance_mat.Prrk1=cell(Nmax,1);
covariance_mat.Prmk1=cell(Nmax,1);
covariance_mat.Pmmk1=cell(Nmax,1);

covariance_mat.Prm_augmented=cell(Nmax,1);
covariance_mat.Pmm_augmented=cell(Nmax,1);

%% Model 3D Parameters:

% Import STL:
Cassini_points=stlread('Envisat.stl');

% rescaling of model to fit inside KOZ:
max_dist=max(vecnorm(Cassini_points.Points,2,2));
scalefact=10; % riscalato per avere 20m max dimension:

% Number of points in model reduction:
fv=triangulation(Cassini_points.ConnectivityList,Cassini_points.Points/max_dist*scalefact);
p=trisurf(fv,'visible','off');
% riduzione del numero di punti del modello per avere uno ridotto da cui
% simulare misure (impossibile stesso livello di accuratezza, 17000 punti
% di cui alcuni sovrapposti, appesantirebbe solo i calcoli:
reducepatch(p,0.2);
fv2=triangulation(p.Faces,p.Vertices);
close all % sennò tiene figure aperta

% Set of Points used by Simulator:
S=fv2.Points;

MASK=zeros(size(fv2.Points,1),size(fv2.ConnectivityList,1));
V = vertexAttachments(fv2);
for i=1:size(V,1)
    MASK(i,V{i})=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cam params:
cam_params.foc=50*1e-3;
cam_params.vpix=1944;
cam_params.hpix=2048;
cam_params.u0=0;
cam_params.v0=0;
cam_params.Hf=0.014566094290851;
cam_params.Vf=cam_params.Hf*cam_params.vpix/cam_params.hpix;
cam_params.b=1;
cam_params.alpha_v=cam_params.foc*cam_params.vpix/cam_params.Vf;
cam_params.alpha_u=cam_params.foc*cam_params.hpix/cam_params.Hf;


%% orbita
params.mu=3.986*1e14;
t0=240*pi/180;
% t0=rand(1);
orbit.RAAN=22*pi/180;
% orbit.i=88.5*pi/180;
orbit.i=98.18*pi/180;
orbit.arg_peri=0*rand(1)*pi/180;
f0=100*pi/180;
a=7138*1e3;
% rp=(6378+700)*1000;         % pericenter
ec=0.0;                    % eccentricity
% a=rp/(1-ec);                % semimajor axis
P=a*(1-ec^2);               % P parameter
r0=P/(1+ec*cos(f0));        % initial orbital sadius
df0=sqrt(P*params.mu)/r0^2;        % initial derivative of true anomaly
dr0=sqrt(params.mu/P)*ec*sin(f0);  % initial radial speed
T=2*pi*sqrt(a^3/params.mu);        % orbital period


%% initialization posizione target rispetto al chaser in LVLH
xp0=-0.002;
yp0=-31.17;
zp0=0.0;
vxp0=-3.5e-6;    % relative speed in x
vyp0=-2.0e-6;       % relative speed in y (determined to be elliptical orbit if e=0)
% vy0=-2.0*x0*(2*pi/T);
vzp0=0.0;    % relative speed in z
s0=[-0.367;-0.590;-0.570];
omega0=[0.02;0.02;0.04];
Jx=125;
Jy=96;
Jz=100;
k10=log(Jx/Jy);
k20=log(Jy/Jz);


%% rotations
d=sqrt(xp0^2+yp0^2+zp0^2);
b0=(-zp0+d)/sqrt(2*d*(-zp0+d));
b1=yp0/sqrt(2*d*(-zp0+d));
b2=-xp0/sqrt(2*d*(-zp0+d));
b3=0;
C_CL=quat2dcm([b0,b1,b2,b3]);
C_LC=C_CL';
C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
qc=dcm2quat(C_CL*C_LI)';
sc0=qc(2:end)/(1+qc(1));
skew_s0=[0,-s0(3),s0(2);s0(3),0,-s0(1);-s0(2),s0(1),0];
% define DCM that brings from chaser to target frame D:
C_TC=eye(3)+8*(skew_s0*skew_s0)/(1+transpose(s0)*s0)^2-4*(1-transpose(s0)*s0)/(1+transpose(s0)*s0)^2*skew_s0;



%% chaser params

Ix=1/12*12*(0.2^2+0.1^2);
Iy=1/12*12*(0.3^2+0.1^2);
Iz=1/12*12*(0.3^2+0.2^2);
I_chas=diag([Ix;Iy;Iz]);
params.Ic_mat=I_chas;

omega_r=cross([xp0;yp0;zp0],[vxp0;vyp0;vzp0])/norm([xp0;yp0;zp0])^2;

omega_r_chas=C_CL*omega_r;
omega_lvlh_i=C_CL*[0;0;df0];
wc0=omega_r_chas+omega_lvlh_i;


x0=[xp0;yp0;zp0;vxp0;vyp0;vzp0;r0;dr0;f0;df0;omega0;s0;k10;k20;wc0;sc0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncertainty Parameters:

% K è moltiplicatore di covarianza, K=0.1,1,3,5,10
K=1;
% Da tesi Pesce
% angular speed:
sr=K*0.1;
sv=K*0.1;
ss=K*0.01;
sw=K*0.01;
sk=K*sqrt(0.1);

% Camera noise, normal distribution
sp=1/sqrt(12);
sp=0.17;

%% Filter Parameters:

% Measure noise covariance:
R=sparse(sp^2*eye(3));

% Process noise covariance:
Q=sparse(1e-7*eye(14));

% Propagation timestep between measures:
dt=0.5;

%% State and Covariance Initialization
% 
% % State Initialization:
% % initial error = 2 sigma:
erk0=1;
err_r=erk0*sr*(2*rand(3,1)-1);
err_v=2*erk0*sv*(2*rand(3,1)-1);
err_w=erk0*sw*(2*rand(3,1)-1);
err_s=erk0*ss*(2*rand(3,1)-1);
err_k=erk0*sk*(2*rand(2,1)-1);
er=[err_r;err_v;zeros(4,1);err_w;err_s;err_k;zeros(6,1)];
spr=0.02;
% % er=zeros(18,1);
% % Initial state is X0+error:

% Initial Point maps:
% S0=fv2.Points(1:5,:)+normrnd(0,spr,5,3);
S0=zeros(0,3);

% Covariance Initialization:
P0=sparse(eye(14+3*size(S0,1)));
P0(1:3,1:3)=sr^2*P0(1:3,1:3);
P0(4:6,4:6)=sv^2*P0(4:6,4:6);
P0(7:9,7:9)=sw^2*P0(7:9,7:9);
P0(10:12,10:12)=ss^2*P0(10:12,10:12);
P0(13:14,13:14)=sk^2*P0(13:14,13:14);

% X0=x0+er;
X0=x0;
% X0([1:6,11:18])=mvnrnd(x0([1:6,11:18]),P0);


Prr0=P0(1:14,1:14);
Prm0=sparse(zeros(14,3*size(S0,1)));
Pmm0=sparse(zeros(3*size(S0,1)));

map_visibility=linspace(1,size(S0,1),size(S0,1))';

%% Plots initialization

figure;
sgtitle('Error evolution')
hAxes(1)=subplot(2,3,1);
hAxes(2)=subplot(2,3,2);
hAxes(4)=subplot(2,3,[3,6]);
hAxes(3)=subplot(2,3,4);
hAxes(5)=subplot(2,3,5);

h1x = animatedline(hAxes(1),'color',[0 0.4470 0.7410]);
h1y = animatedline(hAxes(1),'color',[0.8500 0.3250 0.0980]);
h1z = animatedline(hAxes(1),'color',[0.9290 0.6940 0.1250]);
h1r = animatedline(hAxes(1),'color',[0.4940 0.1840 0.5560]);
h1mdr = animatedline(hAxes(1),'color',[0.4660 0.6740 0.1880]);


%legend:
h1x.DisplayName='dx';
h1y.DisplayName='dy';
h1z.DisplayName='dz';
h1r.DisplayName='dr';
h1mdr.DisplayName='3\sigma_r';
legend(hAxes(1));
grid(hAxes(1),'minor');
xlabel(hAxes(1),'sample n');
ylabel(hAxes(1),'error [m]');


h2x = animatedline(hAxes(2),'color',[0 0.4470 0.7410]);
h2y = animatedline(hAxes(2),'color',[0.8500 0.3250 0.0980]);
h2z = animatedline(hAxes(2),'color',[0.9290 0.6940 0.1250]);
h2r = animatedline(hAxes(2),'color',[0.4940 0.1840 0.5560]);
h2mdr = animatedline(hAxes(2),'color',[0.4660 0.6740 0.1880]);

%legend:
h2x.DisplayName='dvx';
h2y.DisplayName='dvy';
h2z.DisplayName='dvz';
h2r.DisplayName='dvr';
h2mdr.DisplayName='3\sigma_{vr}';
legend(hAxes(2));
grid(hAxes(2),'minor');
xlabel(hAxes(2),'sample n');
ylabel(hAxes(2),'error [m/s]');


h3x = animatedline(hAxes(3),'color',[0 0.4470 0.7410]);
h3y = animatedline(hAxes(3),'color',[0.8500 0.3250 0.0980]);
h3z = animatedline(hAxes(3),'color',[0.9290 0.6940 0.1250]);
h3r = animatedline(hAxes(3),'color',[0.4940 0.1840 0.5560]);
h3mdr = animatedline(hAxes(3),'color',[0.4660 0.6740 0.1880]);

%legend:
h3x.DisplayName='dwx';
h3y.DisplayName='dwy';
h3z.DisplayName='dwz';
h3r.DisplayName='dwr';
h3mdr.DisplayName='3\sigma_{wr}';
legend(hAxes(3));
grid(hAxes(3),'minor');
xlabel(hAxes(3),'sample n');
ylabel(hAxes(3),'error [rad/s]');

h41 = animatedline(hAxes(4),'color',[0 0.4470 0.7410]);
h42 = animatedline(hAxes(4),'color',[0.8500 0.3250 0.0980]);
h4md1 = animatedline(hAxes(4),'color',[0.9290 0.6940 0.1250]);
h4md2 = animatedline(hAxes(4),'color',[0.4940 0.1840 0.5560]);

%legend:
h41.DisplayName='dk1';
h42.DisplayName='dk2';
h4md1.DisplayName='3\sigma_{k1}';
h4md2.DisplayName='3\sigma_{k2}';
legend(hAxes(4));
grid(hAxes(4),'minor');
xlabel(hAxes(4),'sample n');
ylabel(hAxes(4),'error [-]');

h5 = animatedline(hAxes(5),'color',[0 0.4470 0.7410]);

%legend:
h5.DisplayName='de';
legend(hAxes(5));
grid(hAxes(5),'minor');
xlabel(hAxes(5),'sample n');
ylabel(hAxes(5),'error [rad]');