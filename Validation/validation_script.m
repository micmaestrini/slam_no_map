%% Validation Script
% Nmax=159
%% before correction:
% relative position error before correction:
dx_=simulator_history.state(1,:)-state_history.xk(1,:);
dy_=simulator_history.state(2,:)-state_history.xk(2,:);
dz_=simulator_history.state(3,:)-state_history.xk(3,:);

dr=vecnorm([dx_;dy_;dz_],2,1);
pos_ekf=zeros(3,Nmax);
pos_sim=zeros(3,Nmax);
md_r=zeros(Nmax,1);
for i =1:Nmax
cov_r=covariance_mat.Prrk{i}(1:3,1:3);
dr_x=[2*dx_(i);2*dy_(i);2*dz_(i)]'/dr(i);
P_cri=dr_x*cov_r*dr_x';
md_r(i)=3*sqrt(P_cri);

f0=simulator_history.state(9,i);
C_BI=quat2dcm(qc);
C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
C_BL=C_BI*C_LI';
pos_ekf(:,i)=C_BL*[state_history.xk1(1,i);state_history.xk1(2,i);state_history.xk1(3,i)];
pos_sim(:,i)=C_BL*[simulator_history.state(1,i);simulator_history.state(2,i);simulator_history.state(3,i)];
end


% relative velocity error before correction:
dvx_=simulator_history.state(4,:)-state_history.xk(4,:);
dvy_=simulator_history.state(5,:)-state_history.xk(5,:);
dvz_=simulator_history.state(6,:)-state_history.xk(6,:);

dvr=vecnorm([dvx_;dvy_;dvz_],2,1);

md_vr=zeros(Nmax,1);
for i =1:Nmax
cov_vr=covariance_mat.Prrk{i}(4:6,4:6);
dvr_x=[2*dvx_(i);2*dvy_(i);2*dvz_(i)]'/dvr(i);
P_cvri=dvr_x*cov_vr*dvr_x';
md_vr(i)=3*sqrt(P_cvri);
end

% relative angular velocity error before correction:
dwx_=simulator_history.state(11,:)-state_history.xk(11,:);
dwy_=simulator_history.state(12,:)-state_history.xk(12,:);
dwz_=simulator_history.state(13,:)-state_history.xk(13,:);

dwr=vecnorm([dwx_;dwy_;dwz_],2,1);

md_wr=zeros(Nmax,1);
for i =1:Nmax
cov_wr=covariance_mat.Prrk{i}(7:9,7:9);
dwr_x=[2*dwx_(i);2*dwy_(i);2*dwz_(i)]'/dwr(i);
P_cwri=dwr_x*cov_wr*dwr_x';
md_wr(i)=3*sqrt(P_cwri);
end

% relative inertia perameters error before correction:
dk_1=simulator_history.state(17,:)-state_history.xk(17,:);
dk_2=simulator_history.state(18,:)-state_history.xk(18,:);

dk=vecnorm([dk_1;dk_2],2,1);

md_kr=zeros(Nmax,1);
for i =1:Nmax
cov_kr=covariance_mat.Prrk{i}(13:14,13:14);
dkr_x=[2*dk_1(i);2*dk_2(i)]'/dk(i);
P_ckri=dkr_x*cov_kr*dkr_x';
md_kr(i)=3*sqrt(P_ckri);
end

% relative attitude error before correction:
s_sim=simulator_history.state(14:16,:)';
s_prop=state_history.xk(14:16,:)';

ss2=vecnorm(s_sim,2,2).^2;
sp2=vecnorm(s_prop,2,2).^2;

q_sim0=(1-ss2)./(1+ss2);
q_prop0=(1-sp2)./(1+sp2);

q_sim=2*s_sim./(1+ss2);
q_prop=2*s_prop./(1+sp2);

quatsim=[q_sim0,q_sim];
quatprop=[q_prop0,q_prop];

qe=quatmultiply(quatsim,quatinv(quatprop));

err=2*acos(qe(:,1));


%% after correction:

% relative position error before correction:
dx_1=simulator_history.state(1,:)-state_history.xk1(1,:);
dy_1=simulator_history.state(2,:)-state_history.xk1(2,:);
dz_1=simulator_history.state(3,:)-state_history.xk1(3,:);

dr1=vecnorm([dx_1;dy_1;dz_1],2,1);

md_r1=zeros(Nmax,1);
for i =1:Nmax
cov_r1=covariance_mat.Prrk1{i}(1:3,1:3);
dr1_x=[2*dx_1(i);2*dy_1(i);2*dz_1(i)]'/dr1(i);
P_cr1i=dr1_x*cov_r1*dr1_x';
md_r1(i)=3*sqrt(P_cr1i);
end


% relative velocity error before correction:
dvx_1=simulator_history.state(4,:)-state_history.xk1(4,:);
dvy_1=simulator_history.state(5,:)-state_history.xk1(5,:);
dvz_1=simulator_history.state(6,:)-state_history.xk1(6,:);

dvr1=vecnorm([dvx_1;dvy_1;dvz_1],2,1);

md_vr1=zeros(Nmax,1);
for i =1:Nmax
cov_vr1=covariance_mat.Prrk1{i}(4:6,4:6);
dvr1_x=[2*dvx_1(i);2*dvy_1(i);2*dvz_1(i)]'/dvr1(i);
P_cvr1i=dvr1_x*cov_vr1*dvr1_x';
md_vr1(i)=3*sqrt(P_cvr1i);
end

% relative angular velocity error before correction:
dwx_1=simulator_history.state(11,:)-state_history.xk1(11,:);
dwy_1=simulator_history.state(12,:)-state_history.xk1(12,:);
dwz_1=simulator_history.state(13,:)-state_history.xk1(13,:);

dwr1=vecnorm([dwx_1;dwy_1;dwz_1],2,1);

md_wr1=zeros(Nmax,1);
for i =1:Nmax
cov_wr1=covariance_mat.Prrk1{i}(7:9,7:9);
dwr1_x=[2*dwx_1(i);2*dwy_1(i);2*dwz_1(i)]'/dwr1(i);
P_cwr1i=dwr1_x*cov_wr1*dwr1_x';
md_wr1(i)=3*sqrt(P_cwr1i);
end


% relative inertia perameters error before correction:
dk_11=simulator_history.state(17,:)-state_history.xk1(17,:);
dk_21=simulator_history.state(18,:)-state_history.xk1(18,:);

dk1=vecnorm([dk_11;dk_21],2,1);

md_kr1=zeros(Nmax,1);
for i =1:Nmax
cov_kr1=covariance_mat.Prrk1{i}(13:14,13:14);
dkr1_x=[2*dk_11(i);2*dk_21(i)]'/dk1(i);
P_ckr1i=dkr1_x*cov_kr1*dkr1_x';
md_kr1(i)=3*sqrt(P_ckr1i);
end

% relative attitude error before correction:
s_sim1=simulator_history.state(14:16,:)';
s_prop1=state_history.xk1(14:16,:)';

ss21=vecnorm(s_sim1,2,2).^2;
sp21=vecnorm(s_prop1,2,2).^2;

q_sim01=(1-ss21)./(1+ss21);
q_prop01=(1-sp21)./(1+sp21);

q_sim1=2*s_sim1./(1+ss21);
q_prop1=2*s_prop1./(1+sp21);

quatsim1=[q_sim01,q_sim1];
quatprop1=[q_prop01,q_prop1];

qe1=quatmultiply(quatsim1,quatinv(quatprop1));

err1=2*acos(qe1(:,1));
