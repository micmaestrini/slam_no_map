%% parameters initialization tool
% random variables parameter:
clear
clc
close all

%% add paths for functions
addpath('Partials');
addpath('Model');
addpath('Validation');
addpath('Matching_JCBB');

rng(124)

Nmax=300;

times_history=struct('iter',zeros(Nmax,5));

simulator_history=struct('state',zeros(18,Nmax));

simulator_history.measures=cell(Nmax,1);

state_history=struct('xk',zeros(18,Nmax),'xk1',zeros(18,Nmax));

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
scalefact=2; % riscalato per avere 20m max dimension:

% Number of points in model reduction:
fv=triangulation(Cassini_points.ConnectivityList,Cassini_points.Points/max_dist*scalefact);
p=trisurf(fv,'visible','off');
% riduzione del numero di punti del modello per avere uno ridotto da cui
% simulare misure (impossibile stesso livello di accuratezza, 17000 punti
% di cui alcuni sovrapposti, appesantirebbe solo i calcoli:
reducepatch(p,0.002);
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
cam_params.foc=200*1e-3;
cam_params.vpix=400;
cam_params.hpix=600;
cam_params.u0=0;
cam_params.v0=0;
cam_params.Hf=36*1e-3;
cam_params.Vf=24*1e-3;
cam_params.b=1.5;
cam_params.alpha_v=cam_params.foc*cam_params.vpix/cam_params.Vf;
cam_params.alpha_u=cam_params.foc*cam_params.hpix/cam_params.Hf;


%% orbita
mu=3.986*1e14;
t0=10800000;
% t0=rand(1);
orbit.RAAN=30*pi/180;
% orbit.i=88.5*pi/180;
orbit.i=20*pi/180;
orbit.arg_peri=41.4*pi/180;
f0=rand(1)*2*pi;
rp=(6378+700)*1000;         % pericenter
ec=0.02;                    % eccentricity
a=rp/(1-ec);                % semimajor axis
P=a*(1-ec^2);               % P parameter
r0=P/(1+ec*cos(f0));        % initial orbital sadius
df0=sqrt(P*mu)/r0^2;        % initial derivative of true anomaly
dr0=sqrt(mu/P)*ec*sin(f0);  % initial radial speed
T=2*pi*sqrt(a^3/mu);        % orbital period



%% initialization posizione target rispetto al chaser in LVLH
xp0=-100*(2*rand(1)-1);
yp0=-100*(2*rand(1)-1);
zp0=-100*(2*rand(1)-1);
vxp0=0.01*(rand(1)-0.5);    % relative speed in x
vyp0=-2*xp0*(2*pi/T);       % relative speed in y (determined to be elliptical orbit if e=0)
vzp0=0.01*(rand(1)-0.5);    % relative speed in z
s0=rand(3,1);
omega0=rand(3,1)-0.5;
omega0=1*omega0/norm(omega0)*(3*pi/180*rand(1));
Jx=125;
Jy=96;
Jz=100;
k10=log(Jx/Jy);
k20=log(Jy/Jz);
x0=[xp0;yp0;zp0;vxp0;vyp0;vzp0;r0;dr0;f0;df0;omega0;s0;k10;k20];


%% rotations
d=sqrt(xp0^2+yp0^2+zp0^2);
b0=(-zp0+d)/sqrt(2*d*(-zp0+d));
b1=yp0/sqrt(2*d*(-zp0+d));
b2=-xp0/sqrt(2*d*(-zp0+d));
b3=0;
C_CL=quat2dcm([b0,b1,b2,b3]);
C_LC=C_CL';
C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
qc=dcm2quat(C_CL*C_LI);
skew_s0=[0,-s0(3),s0(2);s0(3),0,-s0(1);-s0(2),s0(1),0];
% define DCM that brings from chaser to target frame D:
C_TC=eye(3)+8*(skew_s0*skew_s0)/(1+transpose(s0)*s0)^2-4*(1-transpose(s0)*s0)/(1+transpose(s0)*s0)^2*skew_s0;

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
dt=0.2;

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
er=[err_r;err_v;zeros(4,1);err_w;err_s;err_k];
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

X0=x0+er;
% X0=x0;
% X0([1:6,11:18])=mvnrnd(x0([1:6,11:18]),P0);


Prr0=P0(1:14,1:14);
Prm0=sparse(zeros(14,3*size(S0,1)));
Pmm0=sparse(zeros(3*size(S0,1)));

map_visibility=linspace(1,size(S0,1),size(S0,1))';