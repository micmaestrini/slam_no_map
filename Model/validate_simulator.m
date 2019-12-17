clear
clc
close all

rng(420)

Nmax=400;

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
reducepatch(p,0.1);
fv2=triangulation(p.Faces,p.Vertices);
close all % sennò tiene figure aperta

% Set of Points used by Simulator:
S=fv2.Points;

MASK=zeros(size(fv2.Points,1),size(fv2.ConnectivityList,1));
V = vertexAttachments(fv2);
for i=1:size(V,1)
    MASK(i,V{i})=1;
end

%% orbital Parameters:
% SI units (s, m, kg)

% orbit around Eart:
mu=3.986*1e14;

% Elliptical orbit parameters (tesi Pesce):
a=(7170)*1000;         % pericenter
ec=0.05;                    % eccentricity
% a=rp/(1-ec);                % semimajor axis
P=a*(1-ec^2);               % P parameter
f0=340*pi/180;               % initial anomaly on chaser orbit
r0=P/(1+ec*cos(f0));        % initial orbital sadius
df0=sqrt(P*mu)/r0^2;        % initial derivative of true anomaly
dr0=sqrt(mu/P)*ec*sin(f0);  % initial radial speed
T=2*pi*sqrt(a^3/mu);        % orbital period

% Relative orbit parameters:

xp0=10;      % relative x (aligned with radial)
yp0=60;      % relative y (right hand rule)
zp0=10;      % relative z (aligned with orbital momentum)
vxp0=0.01;    % relative speed in x
vyp0=-0.0225;       % relative speed in y (determined to be elliptical orbit if e=0)
vzp0=-0.01;    % relative speed in z
% MRP of relative attitude:
s0=rand(3,1);

% Relative angular speed C-T (max 3°/s):
% omega0=rand(3,1)-0.5;
% omega0=omega0/norm(omega0)*(3*pi/180*rand(1));
omega0=1*[-0.1;-0.1;0.034]*pi/180;
% Moments of inertia target with parametrization:
Jx=125;
Jy=96;
Jz=100;
k10=log(0.83);
k20=-log(1.083);

% initial state of simulator:
x0=[xp0;yp0;zp0;vxp0;vyp0;vzp0;r0;dr0;f0;df0;omega0;s0;k10;k20];

%% Camera Parameters:

% baseline (1.5 m) along y axis:
cam_params.b=1.0;

% focal length 50mm TBD:
cam_params.foc=30*1e-3;

% Size of sensor (Horizontal and Vertical):
cam_params.Vf=24*1e-3;
cam_params.Hf=36*1e-3;

% Number of pizels in short size:
cam_params.vpix=4000;
cam_params.hpix=4000*1.5;

% pixel densities
cam_params.alpha_v=cam_params.foc*cam_params.vpix/cam_params.Vf;
cam_params.alpha_u=cam_params.foc*cam_params.hpix/cam_params.Hf;

% center pixel offset:
cam_params.u0=0;
cam_params.v0=0;

% initial quaternion (ruota da L a B)
d=sqrt(xp0^2+yp0^2+zp0^2);
b0=(zp0+d)/sqrt(2*d*(zp0+d));
b1=-yp0/sqrt(2*d*(zp0+d));
b2=xp0/sqrt(2*d*(zp0+d));
b3=0;
sc=[b1;b2;b3]/(1+b0);
skew_sc=[0,-sc(3),sc(2);sc(3),0,-sc(1);-sc(2),sc(1),0];
C_BL=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;

C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
qc=dcm2quat(C_BL*C_LI);
        

%% Uncertainty Parameters:

% K è moltiplicatore di covarianza, K=0.1,1,3,5,10
K=1;
% Da tesi Pesce
% angular speed:
sr=K*1;
sv=K*1;
ss=0;
sw=0;
sk=0;

% Camera noise, normal distribution
sp=1/sqrt(12);


%% Filter Parameters:

% Measure noise covariance:
R=sparse(sp^2*eye(3));

% Process noise covariance:
Q=sparse(zeros(6));

% Propagation timestep between measures:
dt=1;

%% State and Covariance Initialization
% 
% % State Initialization:
% % initial error = 2 sigma:
% erk0=2;
% err_r=erk0*sr*(2*rand(3,1)-1);
% err_v=erk0*sv*(2*rand(3,1)-1);
% err_w=erk0*sw*(2*rand(3,1)-1);
% err_s=erk0*ss*(2*rand(3,1)-1);
% err_k=erk0*sk*(2*rand(2,1)-1);
% er=[err_r;err_v;zeros(4,1);err_w;err_s;err_k];
spr=0.02;
% % er=zeros(18,1);
% % Initial state is X0+error:
% X0=x0+er;
% Initial Point maps:
% S0=fv2.Points(1:5,:)+normrnd(0,spr,5,3);
S0=zeros(0,3);

% Covariance Initialization:
P0=sparse(eye(6+3*size(S0,1)));
P0(1:3,1:3)=sr^2*P0(1:3,1:3);
P0(4:6,4:6)=sv^2*P0(4:6,4:6);

X0=x0;
X0(1:6)=mvnrnd(x0(1:6),P0(1:6,1:6));


Prr0=P0(1:6,1:6);
Prm0=sparse(zeros(6,3*size(S0,1)));
Pmm0=sparse(zeros(3*size(S0,1)));

map_visibility=linspace(1,size(S0,1),size(S0,1))';

[yn]=test(X0,fv2,fv2,cam_params,qc,MASK)
% [xn,yn]=simulate_next_step(x0,dt,fv,fv2,cam_params,qc);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % This function simulates a time step in the simulator.
    % Inputs:
    % X    :  state of simulator after dt;
    % fv   :  triangulation object of detailed model;
    % fv2  :  triangulation object of reduced model;
    % npix :  number of pixels on vertical side of sensor;
    % Hf   :  horizontal size of sensor [m];
    % Vf   :  vertical size of sensor [m];
    % foc  :  focal length for camera [m];
    % b    :  baseline between cameras along y [m];
    %
    % Outputs:
    % y    : vector of size [nh,m] of m measurements of size nh each;

    function [y]=test(X,fv,fv2,cam_params,qc,MASK)
    %% conversion of points from target to chaser frame:
    % fv2 for points, fv for triangles
    % retrieve position and relative orientation from current simulator state
    % vector:
        x0=X(1);
        y0=X(2);
        z0=X(3);
        theta=X(9);
        s=X(14:16);
        s1=s(1);
        s2=s(2);
        s3=s(3);

    % define skew symmetric matrix as cross(s,.).
        skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
    % define DCM that brings from chaser to target frame D:
        D=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;
    % D' is DCM from target to chaser:
        D=D';
    % S is the array containing all points in the complete model:
        S=fv.Points;

    % definition of the DCM that rotates from LVLH to chaser body frame:
%     sc=[sc1;sc2;sc3];
%         skew_sc=[0,-sc(3),sc(2);sc(3),0,-sc(1);-sc(2),sc(1),0];
%         Dc=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;
%         C_BL=simplify(Dc);
        C_BI=quat2dcm(qc);
        C_LI=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
        C_BL=C_BI*C_LI';
        qc=dcm2quat(C_BL);
        sc=qc(2:end)/(1+qc(1));
        

    % definition of points of the complete and reduced model in the camera
    % frame (i.e. aligned with x axis of chaser):
        P_i=D*S'+C_BL*[x0;y0;z0];
        P_i2=D*fv2.Points'+C_BL*[x0;y0;z0];

    % definition of the rotated triangulations:
        fv_new=triangulation(fv.ConnectivityList,P_i');
        fv_new2=triangulation(fv2.ConnectivityList,P_i2');

    % extraction of the first, second, and third point of each triangle from
    % the detailed model. extracted as arrays for GPU implementation:
        P0=fv_new2.Points(fv_new2.ConnectivityList(:,1),:);
        P1=fv_new2.Points(fv_new2.ConnectivityList(:,2),:);
        P2=fv_new2.Points(fv_new2.ConnectivityList(:,3),:);

    % definition of matrix of directions from 0 (chaser center of mass) to each
    % of the reduced model points. Also definition of directions from second
    % camera (placed at -b along y axis):
        mat_D=fv_new2.Points;
        mat_D2=mat_D-[cam_params.b,0,0];
        pnorm=vecnorm(mat_D,2,2);
        pnorm2=vecnorm(mat_D2,2,2);
        mat_D=mat_D./pnorm;
        mat_D2=mat_D2./pnorm2;

    % conversion in gpuarrays for gpu implementation:
        gD=gpuArray(mat_D);
        gD2=gpuArray(mat_D2);

    %% Ray tracing step:
    % gpu array function that detects intersection of infinite rays with
    % defined triangles for bot left and right camera:
        [dist, flag] = arrayfun(@rayTriGPU, P0(:,1)', P0(:,2)', P0(:,3)', ...
        P1(:,1)', P1(:,2)', P1(:,3)', ...
        P2(:,1)', P2(:,2)', P2(:,3)', ...
        0, 0, 0, ...
        gD(:,1),gD(:,2),gD(:,3));
        [dist2, flag2] = arrayfun(@rayTriGPU, P0(:,1)', P0(:,2)', P0(:,3)', ...
        P1(:,1)', P1(:,2)', P1(:,3)', ...
        P2(:,1)', P2(:,2)', P2(:,3)', ...
        cam_params.b, 0, 0, ...
        gD2(:,1),gD2(:,2),gD2(:,3));

    %% processing of intersections for both cameras:
    % conversion to normal array (from gpu):
        distances = gather(dist);
    % substitute NaN with inf in array of distances from origin of rays:
        distances(~isfinite(distances)) = inf;
        distances(find(MASK))=Inf;
    % for each ray extract the minimum distance of intersection:
        minint=min(distances,[],2); 

    % conversion to normal array (from gpu):
        distances2 = gather(dist2);
    % substitute NaN with inf in array of distances from origin of rays. Nan
    % values appear if no intersection is found (e.g. point is perfectly
    % visible and is not in front or behind anything).
        distances2(~isfinite(distances2)) = Inf;
        distances2(find(MASK))=Inf;
    % for each ray extract the minimum distance of intersection:
        minint2=min(distances2,[],2);

    %% visibility check:
    % in order for the feature to be visible it must be visible for both
    % cameras. Means that the minimum distance from intersection must be larger
    % than the distance from the camera (i.e. feature is visible but in front
    % of one of the triangular surfaces:
        vis_index=find(pnorm2<minint2 & pnorm<minint);

        points=fv_new2.Points(pnorm2<minint2 & pnorm<minint,:);
        
        trisurf(fv_new2);
        hold on;
        scatter3(fv_new2.Points(:,1),fv_new2.Points(:,2),fv_new2.Points(:,3),'LineWidth',2);
        scatter3(points(:,1),points(:,2),points(:,3),'LineWidth',2);
        campos([0,0,0]);
        
        
        
        
        
        %%%%% WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         points=fv_new2.Points;
%         vis_index=linspace(1,size(points,1),size(points,1))';

    %% measurements conversion to pixels:
        P_i=points';
        xi=P_i(1,:);
        yi=P_i(2,:);
        zi=P_i(3,:);
    % stereo conversion measurement equation:
        h=[cam_params.u0+cam_params.alpha_u*xi./zi;cam_params.v0+cam_params.alpha_v*yi./zi;cam_params.alpha_u*cam_params.b./zi];

%         H_fun(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,fv2.Points(1,1),fv2.Points(1,2),fv2.Points(1,3),s1,s2,s3,sc(1),sc(2),sc(3),cam_params.u0,cam_params.v0,x0,y0,z0)
        %         scatter(h(1,:),h(2,:))
    % points are registered by pixels if they are in FOV of camera (between +
    % and - the size of the sensor):
        vis=abs(h(1,:))<cam_params.hpix/2 & abs(h(2,:))<cam_params.vpix/2;
        in_FOV_index=vis_index(vis);
        vis_meas=h(:,vis>0);
%         vis_meas(1:2,:)=fix(vis_meas(1:2,:));
        
%         y=[vis_meas];
        
        
        y.m=size(vis_meas,2);
        y.z=reshape(vis_meas,[],1);
%         scatter(vis_meas(1,1:10),vis_meas(2,1:10))
        y.vis=in_FOV_index;
    end