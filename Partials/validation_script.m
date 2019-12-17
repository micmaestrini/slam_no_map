clear
clc
close all

rng(692)

rp=(6378+700)*1000;
mu=3.986*1e14;

ec=0.02;
a=rp/(1-ec);
P=a*(1-ec^2);
f0=26*pi/180;
r0=P/(1+ec*cos(f0));
df0=sqrt(P*mu)/r0^2;
dr0=sqrt(mu/P)*ec*sin(f0);
t_grid=2*pi*sqrt(a^3/mu);

xp0=100*(rand(1)-0.5);
yp0=100*(rand(1)-0.5);
zp0=100*(rand(1)-0.5);
vxp0=0.1*(rand(1)-0.5);
vyp0=-2*xp0*(2*pi/t_grid);
vzp0=0.1*(rand(1)-0.5);


% calcolo q0?

s0=0.1*rand(3,1);

Jx=125;
Jy=96;
Jz=100;
k10=log(Jx/Jy);
k20=log(Jy/Jz);

x0=[xp0;yp0;zp0;vxp0;vyp0;vzp0;r0;dr0;f0;df0;rand(3,1)*0.1;s0;k10;k20];



F0=STM(x0);
Phi0=F0-eye(14);
F=@(X) validate_process_jacobian(X,x0(7:10));
% [x,~,~,~,~,~,jacobian] = lsqnonlin(F,X0([1:6,11:end]),[],[],optimset('MaxIter',0));
% Phi01=jacobian;
% spy(Phi01)
multip=0.1;
X=x0([1:6,11:end]);
M=X+multip*eye(14);

Phi1=(F(M(:,1))-F(x0([1:6,11:end])))/multip;
Phi2=(F(M(:,2))-F(x0([1:6,11:end])))/multip;
Phi3=(F(M(:,3))-F(x0([1:6,11:end])))/multip;
Phi4=(F(M(:,4))-F(x0([1:6,11:end])))/multip;
Phi5=(F(M(:,5))-F(x0([1:6,11:end])))/multip;
Phi6=(F(M(:,6))-F(x0([1:6,11:end])))/multip;
Phi7=(F(M(:,7))-F(x0([1:6,11:end])))/multip;
Phi8=(F(M(:,8))-F(x0([1:6,11:end])))/multip;
Phi9=(F(M(:,9))-F(x0([1:6,11:end])))/multip;
Phi10=(F(M(:,10))-F(x0([1:6,11:end])))/multip;
Phi11=(F(M(:,11))-F(x0([1:6,11:end])))/multip;
Phi12=(F(M(:,12))-F(x0([1:6,11:end])))/multip;
Phi13=(F(M(:,13))-F(x0([1:6,11:end])))/multip;
Phi14=(F(M(:,14))-F(x0([1:6,11:end])))/multip;

Phi1=[Phi1,Phi2,Phi3,Phi4,Phi5,Phi6,Phi7,Phi8,Phi9,Phi10,Phi11,Phi12,Phi13,Phi14];

Phi0-Phi1
% Phi0(1:6,1:6)-[Phi1(1:6),Phi2(1:6),Phi3(1:6),Phi4(1:6),Phi5(1:6),Phi6(1:6)]
% Phi0(10:13,:)-[Phi1(10:13),Phi2(10:13),Phi3(10:13),Phi4(10:13),Phi5(10:13),Phi6(10:13),Phi7(10:13),Phi8(10:13),Phi9(10:13),Phi10(10:13),Phi11(10:13),Phi12(10:13),Phi13(10:13),Phi14(10:13),Phi15(10:13)]
% Phi0(7:9,:)-[Phi1(7:9),Phi2(7:9),Phi3(7:9),Phi4(7:9),Phi5(7:9),Phi6(7:9),Phi7(7:9),Phi8(7:9),Phi9(7:9),Phi10(7:9),Phi11(7:9),Phi12(7:9),Phi13(7:9),Phi14(7:9),Phi15(7:9)]
%
% Phi0(7:9,7:9)-[Phi7(7:9),Phi8(7:9),Phi9(7:9)]
% Phi0(7:9,10:13)-[Phi10(7:9),Phi11(7:9),Phi12(7:9),Phi13(7:9)]
% Phi0(7:9,14:15)-[Phi14(7:9),Phi15(7:9)]
% Phi0(7:9,14:end)-[Phi14(7:9),Phi15(7:9)]



p=rand(3,1)*1.4;
b=1.5;
foc=200*1e-3;
delta=1e-7;
measures = meas_fun(b,foc,p(1),p(2),p(3),s0(1),s0(2),s0(3),xp0,yp0,zp0);
gm1 = (meas_fun(b,foc,p(1)+delta,p(2),p(3),s0(1),s0(2),s0(3),xp0,yp0,zp0)-measures)/delta;
gm2 = (meas_fun(b,foc,p(1),p(2)+delta,p(3),s0(1),s0(2),s0(3),xp0,yp0,zp0)-measures)/delta;
gm3 = (meas_fun(b,foc,p(1),p(2),p(3)+delta,s0(1),s0(2),s0(3),xp0,yp0,zp0)-measures)/delta;
gm4 = (meas_fun(b,foc,p(1),p(2),p(3),s0(1)+delta,s0(2),s0(3),xp0,yp0,zp0)-measures)/delta;
gm5 = (meas_fun(b,foc,p(1),p(2),p(3),s0(1),s0(2)+delta,s0(3),xp0,yp0,zp0)-measures)/delta;
gm6 = (meas_fun(b,foc,p(1),p(2),p(3),s0(1),s0(2),s0(3)+delta,xp0,yp0,zp0)-measures)/delta;
gm7 = (meas_fun(b,foc,p(1),p(2),p(3),s0(1),s0(2),s0(3),xp0+delta,yp0,zp0)-measures)/delta;
gm8 = (meas_fun(b,foc,p(1),p(2),p(3),s0(1),s0(2),s0(3),xp0,yp0+delta,zp0)-measures)/delta;
gm9 = (meas_fun(b,foc,p(1),p(2),p(3),s0(1),s0(2),s0(3),xp0,yp0,zp0+delta)-measures)/delta;


dh = H(b,foc,p(1),p(2),p(3),s0(1),s0(2),s0(3),xp0,yp0,zp0);
dh_pi = H_pi(b,foc,p(1),p(2),p(3),s0(1),s0(2),s0(3),xp0,yp0,zp0);

dh_app=[gm7,gm8,gm9,zeros(4,6),gm4,gm5,gm6];
dh_pi_app=[gm1,gm2,gm3];

dh-dh_app
dh_pi_app-dh_pi



% (b,foc,s1,s2,s3,x,xl,xr,y,yr,z)
% 
% (b,foc,s1,s2,s3,x,xl,xr,y,yr,z)
% (b,foc,s1,s2,s3,x,xl,xr,y,yr,z)

delta=1e-7;
xr=measures(1);
yr=measures(2);
xl=measures(3);
yl=measures(4);
invmeasures = inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0,xl,xr,yp0,yr,zp0);
gm1 = (inv_meas_fun(b,foc,s0(1)+delta,s0(2),s0(3),xp0,xl,xr,yp0,yr,zp0)-invmeasures)/delta;
gm2 = (inv_meas_fun(b,foc,s0(1),s0(2)+delta,s0(3),xp0,xl,xr,yp0,yr,zp0)-invmeasures)/delta;
gm3 = (inv_meas_fun(b,foc,s0(1),s0(2),s0(3)+delta,xp0,xl,xr,yp0,yr,zp0)-invmeasures)/delta;
gm4 = (inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0+delta,xl,xr,yp0,yr,zp0)-invmeasures)/delta;
gm5 = (inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0,xl+delta,xr,yp0,yr,zp0)-inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0,xl-delta,xr,yp0,yr,zp0))/(2*delta);
gm6 = (inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0,xl,xr+delta,yp0,yr,zp0)-inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0,xl,xr-delta,yp0,yr,zp0))/(2*delta);
gm7 = (inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0,xl,xr,yp0+delta,yr,zp0)-invmeasures)/delta;
gm8 = (inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0,xl,xr,yp0,yr+delta,zp0)-invmeasures)/delta;
gm9 = (inv_meas_fun(b,foc,s0(1),s0(2),s0(3),xp0,xl,xr,yp0,yr,zp0+delta)-invmeasures)/delta;


dg  =    G(b,foc,s0(1),s0(2),s0(3),xp0,xl,xr,yp0,yr,zp0);
dg_pi=G_pi(b,foc,s0(1),s0(2),s0(3),xl,xr,yr);

dg_app=[gm4,gm7,gm9,zeros(3,6),gm1,gm2,gm3];
dg_pi_app=[gm6,gm8,gm5,zeros(3,1)];

dg-dg_app
dg_pi_app-dg_pi

dg_pi(:,1)-gm6
dg_pi(:,3)-gm5









function [dX]=validate_process_jacobian(X,params)
r=params(1);
dr=params(2);
f=params(3);
df=params(4);
mu=3.986*1e14;

x=X(1);
y=X(2);
z=X(3);
dx=X(4);
dy=X(5);
dz=X(6);
w=X(7:9);
s=X(10:12);
k1=X(13);
k2=X(14);

% useful quantities definition
s1=s(1);
s2=s(2);
s3=s(3);
skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
D=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;


% D=quat2dcm(q');


J=diag([exp(k1),1,exp(k2)]);
invJ=diag([1/exp(k1),1,1/exp(k2)]);

B=1/4*((1-transpose(s)*s)*eye(3)+2*skew_s+2*s*transpose(s));

% NERM derivatives of relative orbit
ddf=-2*dr.*df./r;
ddx=2*df.*dy+ddf.*y+df.^2.*x-mu*(r+x)./((r+x).^2+y.^2+z.^2).^1.5+mu./r.^2;
ddy=-2*df.*dx-ddf.*x+df.^2.*y-mu*y./((r+x).^2+y.^2+z.^2).^1.5;
ddz=-mu*z./((r+x).^2+y.^2+z.^2).^1.5;
ddr=r.*df.^2-mu./r.^2;

% derivatives of relative attitude

b = b_fun(x,y,z);
b0=b(1);
b1=b(2);
b2=b(3);
b3=b(4);
db = db_fun(dx,dy,dz,x,y,z);
db0=db(1);
db1=db(2);
db2=db(3);
db3=db(4);
ddb = ddb_fun(df,dr,mu,r,dx,dy,dz,x,y,z);
ddb0=ddb(1);
ddb1=ddb(2);
ddb2=ddb(3);
ddb3=ddb(4);

wc=wc_fun(b0,b1,b2,b3,db0,db1,db2,db3,f,df);
dwc=dwc_fun(b0,b1,b2,b3,db0,db1,db2,db3,ddb0,ddb1,ddb2,ddb3,ddf,f,df);

ds=B*w;
dw=invJ*(-cross(w,J*w)+J*cross(w,D*wc)-J*D*dwc-cross(D*wc,J*D*wc)-cross(w,J*D*wc)-cross(D*wc,J*w));

% derivatives of inertia parametrization
dk1=0;
dk2=0;


dX=[dx;dy;dz;ddx;ddy;ddz;dw;ds;dk1;dk2];

end

