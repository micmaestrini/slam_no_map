function [hn]=new_meas(meas0,X0,Xn,cam_params,indexPairs)

        pix_coord0=reshape(meas0.z,3,[])';
        pix_coord0=pix_coord0(indexPairs(:,1),:);
    % inverse measure function:
        f0=X0(9);
        x0=X0(1);
        y0=X0(2);
        z0=X0(3);
        s0=X0(14:16);
        sc0=X0(22:24);
        
        skew_sc0=[0,-sc0(3),sc0(2);sc0(3),0,-sc0(1);-sc0(2),sc0(1),0];
        C_BI0=eye(3)+8*(skew_sc0*skew_sc0)/(1+transpose(sc0)*sc0)^2-4*(1-transpose(sc0)*sc0)/(1+transpose(sc0)*sc0)^2*skew_sc0;
        C_LI0=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL0=C_BI0*C_LI0';
        
               
        p10=-cam_params.b./pix_coord0(:,3).*(pix_coord0(:,1)-cam_params.u0);
        p20=-cam_params.alpha_u/cam_params.alpha_v*cam_params.b./pix_coord0(:,3).*(pix_coord0(:,2)-cam_params.v0);
        p30=cam_params.alpha_u*cam_params.b./pix_coord0(:,3);
        
        skew_s0=[0,-s0(3),s0(2);s0(3),0,-s0(1);-s0(2),s0(1),0];
        C_TB0=eye(3)+8*(skew_s0*skew_s0)/(1+transpose(s0)*s0)^2-4*(1-transpose(s0)*s0)/(1+transpose(s0)*s0)^2*skew_s0;
        S0=[p10,p20,p30];
        
        P_i0=C_TB0*(S0'-C_BL0*[x0;y0;z0]);
        
        % points in target reference frame at time0:
        Pi0=P_i0';
        
        
        xn=Xn(1);
        yn=Xn(2);
        zn=Xn(3);
        fn=Xn(9);
        sn=Xn(14:16);


    % DCM that rotates from LVLH to chaser body frame definition:
        scn=Xn(22:24);
        skew_scn=[0,-scn(3),scn(2);scn(3),0,-scn(1);-scn(2),scn(1),0];
        C_BIn=eye(3)+8*(skew_scn*skew_scn)/(1+transpose(scn)*scn)^2-4*(1-transpose(scn)*scn)/(1+transpose(scn)*scn)^2*skew_scn;
        C_LIn=[cos(fn),sin(fn),0;-sin(fn),cos(fn),0;0,0,1];
        C_BLn=C_BIn*C_LIn';



%     % DCM that rotates from LVLH to chaser body frame definition:
%         C_BI=quat2dcm(qc);
%         C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
%         C_BL=C_BI*C_LI';
%         qc=dcm2quat(C_BL);
%         sc=qc(2:end)/(1+qc(1));

    % DCM that rotates from target to chaser frame (D'):
        skew_sn=[0,-sn(3),sn(2);sn(3),0,-sn(1);-sn(2),sn(1),0];
        C_TBn=eye(3)+8*(skew_sn*skew_sn)/(1+transpose(sn)*sn)^2-4*(1-transpose(sn)*sn)/(1+transpose(sn)*sn)^2*skew_sn;
        C_BTn=C_TBn';

    % rotation of current landmarks in the chaser reference frame:
        P_in=C_BTn*Pi0'+C_BLn*[xn;yn;zn];
        xin=P_in(1,:);
        yin=P_in(2,:);
        zin=P_in(3,:);

    % stereo measurement equation:
        hn=[cam_params.u0-cam_params.alpha_u*xin./zin;cam_params.v0-cam_params.alpha_v*yin./zin;cam_params.alpha_u*cam_params.b./zin]';
