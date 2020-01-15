clear
clc
close all

%%
        syms x0 y0 z0 f0 s01 s02 s03 sc01 sc02 sc03 u0 v0 d0 b u v alpha_u alpha_v
        
        % first function of the states:
        p10=-b/d0*(u0-u);
        p20=-alpha_u/alpha_v*b/d0*(v0-v);
        p30=alpha_u*b/d0;        
        dm0_dy0=simplify(jacobian([p10;p20;p30],[u0;v0;d0]));
        
        syms m01 m02 m03
        m0=[m01;m02;m03];

        
        %%
    % inverse measure function:
        % not fun of the states:
        sc0=[sc01;sc02;sc03];
        skew_sc0=[0,-sc0(3),sc0(2);sc0(3),0,-sc0(1);-sc0(2),sc0(1),0];
        C_BI0=eye(3)+8*(skew_sc0*skew_sc0)/(1+transpose(sc0)*sc0)^2-4*(1-transpose(sc0)*sc0)/(1+transpose(sc0)*sc0)^2*skew_sc0;
        C_LI0=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL0=simplify(C_BI0*transpose(C_LI0));
            
        % fun of the state:
        s0=[s01;s02;s03];
        skew_s0=[0,-s0(3),s0(2);s0(3),0,-s0(1);-s0(2),s0(1),0];
        C_TB0=eye(3)+8*(skew_s0*skew_s0)/(1+transpose(s0)*s0)^2-4*(1-transpose(s0)*s0)/(1+transpose(s0)*s0)^2*skew_s0;
        C_TL0=simplify(C_TB0*C_BL0);
        
        P_i0=simplify(C_TB0*m0-C_TL0*[x0;y0;z0]);
        
        syms wx0 wy0 wz0 k10 k20 vx0 vy0 vz0
        dP0i_dx0=simplify(jacobian(P_i0,[x0;y0;z0;vx0;vy0;vz0;wx0;wy0;wz0;s0;k10;k20]));
        dP0i_dm0=simplify(jacobian(P_i0,m0));
        
        syms X0i Y0i Z0i
        P0i=[X0i;Y0i;Z0i];
        
        %%        
        syms xn yn zn fn sn1 sn2 sn3 scn1 scn2 scn3

        % not a fun of the states:
        scn=[scn1;scn2;scn3];
        skew_scn=[0,-scn(3),scn(2);scn(3),0,-scn(1);-scn(2),scn(1),0];
        C_BIn=eye(3)+8*(skew_scn*skew_scn)/(1+transpose(scn)*scn)^2-4*(1-transpose(scn)*scn)/(1+transpose(scn)*scn)^2*skew_scn;
        C_LIn=[cos(fn),sin(fn),0;-sin(fn),cos(fn),0;0,0,1];
        C_BLn=simplify(C_BIn*transpose(C_LIn));

    % fun of the states:
        sn=[sn1;sn2;sn3];
        skew_sn=[0,-sn(3),sn(2);sn(3),0,-sn(1);-sn(2),sn(1),0];
        C_TBn=eye(3)+8*(skew_sn*skew_sn)/(1+transpose(sn)*sn)^2-4*(1-transpose(sn)*sn)/(1+transpose(sn)*sn)^2*skew_sn;
        C_BTn=simplify(transpose(C_TBn));
        
    % rotation of current landmarks in the chaser reference frame:
        P_in=simplify(C_BTn*P0i+C_BLn*[xn;yn;zn]);
        
        
        syms wxn wyn wzn k1n k2n vxn vyn vzn
        
        dPni_dxn=simplify(jacobian(P_in,[xn;yn;zn;vxn;vyn;vzn;wxn;wyn;wzn;sn;k1n;k2n]));
        dPni_dP0i=simplify(jacobian(P_in,P0i));
        
%%
    % stereo measurement equation:
    
        syms Xin Yin Zin
        
        hn=[u-alpha_u*Xin/Zin;v-alpha_v*Yin/Zin;alpha_u*b/Zin];
        dhn_dxin=simplify(jacobian(hn,[Xin;Yin;Zin]));
        
        
%% chain rule combinations:

dh_dxn=simplify(dhn_dxin*dPni_dxn);

dh_dxin=dhn_dxin;
dxin_dxi0=dPni_dP0i;
dxi0_dx0=dP0i_dx0;

dh_dy0=simplify(dhn_dxin*dPni_dP0i*dP0i_dm0*dm0_dy0);



%%


matlabFunction(dh_dxn,'File','jacobian_meas_n_state_n','Sparse',true);

matlabFunction(dh_dxin,'File','jacobian_meas_n_points_n','Sparse',true);
matlabFunction(dxin_dxi0,'File','jacobian_points_n_points_0','Sparse',true);
matlabFunction(dxi0_dx0,'File','jacobian_points_0_state_0','Sparse',true);

matlabFunction(dh_dy0,'File','jacobian_meas_n_meas_0','Sparse',true);




        
