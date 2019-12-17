clear
clc
close all

%% NB consider camera left centered in CDM of chaser and aligned with zc axis.
% moreover C is at +b distance along xc axis:

syms sc1 sc2 sc3 x y z vx vy vz wx wy wz
sc=[sc1;sc2;sc3];
skew_sc=[0,-sc3,sc2;sc3,0,-sc1;-sc2,sc1,0];
Dc=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;
C_BL=simplify(Dc);
%%
syms s1 s2 s3 p1 p2 p3 alpha_u alpha_v b k1 k2 u0 v0
s=[s1;s2;s3];
skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
D=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;
D=simplify(D);
D=transpose(D);

P_i=simplify((D*[p1;p2;p3]+C_BL*[x;y;z]));

xi=P_i(1);
yi=P_i(2);
zi=P_i(3);
%%
h=simplify([u0-alpha_u*xi/zi;v0-alpha_v*yi/zi;alpha_u*b/zi]);

dh=(jacobian(h,[x,y,z,vx,vy,vz,wx,wy,wz,s1,s2,s3,k1,k2]));
dh_pi=(jacobian(h,[p1,p2,p3]));
%%
matlabFunction(h,'File','H_fun','Sparse',true);
matlabFunction(dh,'File','H_x','Sparse',true);
matlabFunction(dh_pi,'File','H_Pi','Sparse',true);


