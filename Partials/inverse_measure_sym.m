clear
clc
close all


syms u v u0 v0 alpha_u alpha_v b d


p1=-b/d*(u-u0);
p2=-alpha_u/alpha_v*b/d*(v-v0);
p3=alpha_u*b/d;


S_new=[p1;p2;p3];


syms sc1 sc2 sc3 x y z vx vy vz

sc=[sc1;sc2;sc3];
skew_sc=[0,-sc3,sc2;sc3,0,-sc1;-sc2,sc1,0];
Dc=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;
C_BL=simplify(Dc);


syms s1 s2 s3 wx wy wz vx vy vz k1 k2

s=[s1;s2;s3];
skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
Dt=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;


P_i=simplify(Dt*(S_new-C_BL*[x;y;z]));

dg=simplify(jacobian(P_i,[x,y,z,vx,vy,vz,wx,wy,wz,s1,s2,s3,k1,k2]));
dg_pi=simplify(jacobian(P_i,[u,v,d]));

matlabFunction(P_i,'File','G_fun','Sparse',true);
matlabFunction(dg,'File','G_x','Sparse',true);
matlabFunction(dg_pi,'File','G_Pi','Sparse',true);
