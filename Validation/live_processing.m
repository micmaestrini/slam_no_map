%% Validation Script
% Nmax=159
%% before correction:
% relative position error before correction:
delta_x=xn-Xn1;
dx_=delta_x(1);
dy_=delta_x(2);
dz_=delta_x(3);
dr=norm(delta_x(1:3));
        
cov_r=Prr1(1:3,1:3);
dr_x=[2*dx_;2*dy_;2*dz_]'/dr;
P_cri=dr_x*cov_r*dr_x';
md_r=3*sqrt(P_cri);

%%
% relative velocity error before correction:
dvx_=delta_x(4);
dvy_=delta_x(5);
dvz_=delta_x(6);
dvr=norm(delta_x(4:6));
        
cov_vr=Prr1(4:6,4:6);
dvr_x=[2*dvx_;2*dvy_;2*dvz_]'/dvr;
P_cvri=dvr_x*cov_vr*dvr_x';
md_vr=3*sqrt(P_cvri);


%%
% relative angular velocity error before correction:
dwx_=delta_x(11);
dwy_=delta_x(12);
dwz_=delta_x(13);
dwr=norm(delta_x(11:13));
        
cov_wr=Prr1(7:9,7:9);
dwr_x=[2*dwx_;2*dwy_;2*dwz_]'/dwr;
P_cwri=dwr_x*cov_wr*dwr_x';
md_wr=3*sqrt(P_cwri);

%%
% inertia parameters error before correction:

dk1=delta_x(17);
dk2=delta_x(18);    
md_k1=full(3*sqrt(Prr1(13,13)));
md_k2=full(3*sqrt(Prr1(14,14)));


%%
% relative attitude error before correction:

s_sim=xn(14:16);
s_prop=Xn1(14:16);

% s_sim=[(1-norm(mrps_sim)^2)/(1+norm(mrps_sim)^2);2*mrps_sim/(1+norm(mrps_sim)^2)];
% s_prop=[(1-norm(mrps_prop)^2)/(1+norm(mrps_prop)^2);2*mrps_prop/(1+norm(mrps_prop)^2)];

ss2=norm(s_sim).^2;
sp2=norm(s_prop).^2;

q_sim0=(1-ss2)./(1+ss2);
q_prop0=(1-sp2)./(1+sp2);

q_sim=2*s_sim./(1+ss2);
q_prop=2*s_prop./(1+sp2);

quatsim=[q_sim0;q_sim];
quatprop=[q_prop0;q_prop];

qe=quatmultiply(quatsim',quatinv(quatprop'));

err=2*acos(qe(1));

addpoints(h1x,loop,dx_);
addpoints(h1y,loop,dy_);
addpoints(h1z,loop,dz_);
addpoints(h1r,loop,dr);
addpoints(h1mdr,loop,md_r);

addpoints(h2x,loop,dvx_);
addpoints(h2y,loop,dvy_);
addpoints(h2z,loop,dvz_);
addpoints(h2r,loop,dvr);
addpoints(h2mdr,loop,md_vr);

addpoints(h3x,loop,dwx_);
addpoints(h3y,loop,dwy_);
addpoints(h3z,loop,dwz_);
addpoints(h3r,loop,dwr);
addpoints(h3mdr,loop,md_wr);

addpoints(h41,loop,dk1);
addpoints(h42,loop,dk2);
addpoints(h4md1,loop,md_k1);
addpoints(h4md2,loop,md_k2);

addpoints(h5,loop,err);
drawnow








