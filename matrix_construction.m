function [A,B,C,h_x]=matrix_construction(X0,Xn,Y0,cam_params)
%% data extraction
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn=Xn(1);
yn=Xn(2);
zn=Xn(3);
x0=X0(1);
y0=X0(2);
z0=X0(3);

fn=Xn(9);
f0=X0(9);

sn1=Xn(14);
sn2=Xn(15);
sn3=Xn(16);
s01=X0(14);
s02=X0(15);
s03=X0(16);

scn1=Xn(22);
scn2=Xn(23);
scn3=Xn(24);
sc01=X0(22);
sc02=X0(23);
sc03=X0(24);


A_size=3*size(Y0,1);
B_size=3*size(Y0,1);
C_size=3*size(Y0,1);

A=sparse(zeros(A_size,14));
B=sparse(zeros(B_size,14));
C=sparse(zeros(C_size,C_size));
h_x=zeros(3*size(Y0,1),1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% one shot computations:
% inverse measure function:
% not fun of the states:
sc0=[sc01;sc02;sc03];
skew_sc0=[0,-sc0(3),sc0(2);sc0(3),0,-sc0(1);-sc0(2),sc0(1),0];
C_BI0=eye(3)+8*(skew_sc0*skew_sc0)/(1+transpose(sc0)*sc0)^2-4*(1-transpose(sc0)*sc0)/(1+transpose(sc0)*sc0)^2*skew_sc0;
C_LI0=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
C_BL0=(C_BI0*transpose(C_LI0));

scn=[scn1;scn2;scn3];
skew_scn=[0,-scn(3),scn(2);scn(3),0,-scn(1);-scn(2),scn(1),0];
C_BIn=eye(3)+8*(skew_scn*skew_scn)/(1+transpose(scn)*scn)^2-4*(1-transpose(scn)*scn)/(1+transpose(scn)*scn)^2*skew_scn;
C_LIn=[cos(fn),sin(fn),0;-sin(fn),cos(fn),0;0,0,1];
C_BLn=(C_BIn*transpose(C_LIn));

% fun of the state:
s0=[s01;s02;s03];
skew_s0=[0,-s0(3),s0(2);s0(3),0,-s0(1);-s0(2),s0(1),0];
C_TB0=eye(3)+8*(skew_s0*skew_s0)/(1+transpose(s0)*s0)^2-4*(1-transpose(s0)*s0)/(1+transpose(s0)*s0)^2*skew_s0;
C_TL0=(C_TB0*C_BL0);

sn=[sn1;sn2;sn3];
skew_sn=[0,-sn(3),sn(2);sn(3),0,-sn(1);-sn(2),sn(1),0];
C_TBn=eye(3)+8*(skew_sn*skew_sn)/(1+transpose(sn)*sn)^2-4*(1-transpose(sn)*sn)/(1+transpose(sn)*sn)^2*skew_sn;
C_BTn=(transpose(C_TBn));

x_t_c_T=C_TL0*[x0;y0;z0];
x_t_c_C=C_BLn*[xn;yn;zn];

dxin_dxi0 = jacobian_points_n_points_0(sn1,sn2,sn3);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% things to put in loop:

for i=1:size(Y0,1)
    
    u0=Y0(i,1);
    v0=Y0(i,2);
    d0=Y0(i,3);
    
    
    
    p10=-cam_params.b/d0*(u0-cam_params.u0);
    p20=-cam_params.alpha_u/cam_params.alpha_v*cam_params.b/d0*(v0-cam_params.v0);
    p30=cam_params.alpha_u*cam_params.b/d0;
    P_i0=(C_TB0*[p10;p20;p30]-x_t_c_T);
    P_in=(C_BTn*P_i0+x_t_c_C);
    hn=[cam_params.u0-cam_params.alpha_u*P_in(1)/P_in(3);cam_params.v0-cam_params.alpha_v*P_in(2)/P_in(3);cam_params.alpha_u*cam_params.b/P_in(3)];
    
    dhi_dxn = jacobian_meas_n_state_n(P_i0(1),P_in(1),P_i0(2),P_in(2),P_i0(3),P_in(3),cam_params.alpha_u,cam_params.alpha_v,cam_params.b,fn,scn1,scn2,scn3,sn1,sn2,sn3);
    dh_dy0 = jacobian_meas_n_meas_0(P_in(1),P_in(2),P_in(3),cam_params.alpha_u,cam_params.alpha_v,cam_params.b,d0,s01,s02,s03,sn1,sn2,sn3,cam_params.u0,u0,cam_params.v0,v0);
    
    
    dh_dxin = jacobian_meas_n_points_n(P_in(1),P_in(2),P_in(3),cam_params.alpha_u,cam_params.alpha_v,cam_params.b);
    dxi0_dx0 = jacobian_points_0_state_0(f0,p10,p20,p30,s01,s02,s03,sc01,sc02,sc03,x0,y0,z0);
    
    
    A(1+(i-1)*3:(i)*3,:)=dhi_dxn;
    B(1+(i-1)*3:(i)*3,:)=dh_dxin*dxin_dxi0*dxi0_dx0;
    C(1+(i-1)*3:(i)*3,1+(i-1)*3:(i)*3)=dh_dy0;
    
    h_x(1+(i-1)*3:(i)*3)=hn;
    
end



end