clear
clc
close all


        syms x0 y0 z0 f0 s01 s02 s03 sc01 sc02 sc03 u0 v0 d0 b u v alpha_u alpha_v
        
    % inverse measure function:

        s0=[s01;s02;s03];
        sc0=[sc01;sc02;sc03];
        
        skew_sc0=[0,-sc0(3),sc0(2);sc0(3),0,-sc0(1);-sc0(2),sc0(1),0];
        C_BI0=eye(3)+8*(skew_sc0*skew_sc0)/(1+transpose(sc0)*sc0)^2-4*(1-transpose(sc0)*sc0)/(1+transpose(sc0)*sc0)^2*skew_sc0;
        C_LI0=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL0=C_BI0*transpose(C_LI0);
        
               
        p10=-b/d0*(u0-u);
        p20=-alpha_u/alpha_v*b/d0*(v0-v);
        p30=alpha_u*b/d0;
        
        skew_s0=[0,-s0(3),s0(2);s0(3),0,-s0(1);-s0(2),s0(1),0];
        C_TB0=eye(3)+8*(skew_s0*skew_s0)/(1+transpose(s0)*s0)^2-4*(1-transpose(s0)*s0)/(1+transpose(s0)*s0)^2*skew_s0;
        S0=[p10;p20;p30];
        
        P_i0=simplify(C_TB0*(S0-C_BL0*[x0;y0;z0]));
        
        % points in target reference frame at time0:
        Pi0=P_i0;
        
        syms xn yn zn fn sn1 sn2 sn3 scn1 scn2 scn3

        sn=[sn1;sn2;sn3];
        scn=[scn1;scn2;scn3];

        skew_scn=[0,-scn(3),scn(2);scn(3),0,-scn(1);-scn(2),scn(1),0];
        C_BIn=eye(3)+8*(skew_scn*skew_scn)/(1+transpose(scn)*scn)^2-4*(1-transpose(scn)*scn)/(1+transpose(scn)*scn)^2*skew_scn;
        C_LIn=[cos(fn),sin(fn),0;-sin(fn),cos(fn),0;0,0,1];
        C_BLn=C_BIn*transpose(C_LIn);

    % DCM that rotates from target to chaser frame (D'):
        skew_sn=[0,-sn(3),sn(2);sn(3),0,-sn(1);-sn(2),sn(1),0];
        C_TBn=eye(3)+8*(skew_sn*skew_sn)/(1+transpose(sn)*sn)^2-4*(1-transpose(sn)*sn)/(1+transpose(sn)*sn)^2*skew_sn;
        C_BTn=transpose(C_TBn);

    % rotation of current landmarks in the chaser reference frame:
        P_in=C_BTn*Pi0+C_BLn*[xn;yn;zn];
        xin=P_in(1);
        yin=P_in(2);
        zin=P_in(3);

    % stereo measurement equation:
        hn=[u-alpha_u*xin/zin;v-alpha_v*yin/zin;alpha_u*b/zin];
        
        syms wxn wyn wzn k1n k2n
        
        H=jacobian(hn,[xn;yn;zn;wxn;wyn;wzn;sn1;sn2;sn3;k1n;k2n]);
        
        
