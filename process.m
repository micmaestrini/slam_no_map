function [dX]=process(t,X,params)

mu=params.mu;


x=X(1);
y=X(2);
z=X(3);
dx=X(4);
dy=X(5);
dz=X(6);
r=X(7);
dr=X(8);
f=X(9);
df=X(10);
w=X(11:13);
s=X(14:16);
k1=X(17);
k2=X(18);

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

wc=X(19:21);
dwc=-inv(params.Ic_mat)*(cross(wc,params.Ic_mat*wc));
sc=X(22:24);
skew_sc=[0,-sc(3),sc(2);sc(3),0,-sc(1);-sc(2),sc(1),0];
% Dc=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;
Bc=1/4*((1-transpose(sc)*sc)*eye(3)+2*skew_sc+2*sc*transpose(sc));
dsc=Bc*wc;

ds=B*w;
dw=invJ*(-cross(w,J*w)+J*cross(w,D*wc)-J*D*dwc-cross(D*wc,J*D*wc)-cross(w,J*D*wc)-cross(D*wc,J*w));

% derivatives of inertia parametrization
dk1=0;
dk2=0;


dX=[dx;dy;dz;ddx;ddy;ddz;dr;ddr;df;ddf;dw;ds;dk1;dk2;dwc;dsc];
end
