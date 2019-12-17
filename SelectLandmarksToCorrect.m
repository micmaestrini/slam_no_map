function [to_correct]=SelectLandmarksToCorrect(X,Prr,Prm,Pmm,vis,qc,cam_params)

        meas_uncert=zeros(size(vis,1),1);
        f0=Xn(9);
        s0=Xn(14:16);
        
        C_BI=quat2dcm(qc);
        C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL=C_BI*C_LI';
        qc=dcm2quat(C_BL);
        
        sc=qc(2:end)/(1+qc(1));

        x0=X(1);
        y0=X(2);
        z0=X(3);
        
        ns=3;
        
        for i=1:size(vis,1)
        % partial of the measure wrt the state:
        h_xi=H_x(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,S0(vis(i),1),S0(vis(i),2),S0(vis(i),3),s0(1),s0(2),s0(3),sc(1),sc(2),sc(3),x0,y0,z0);
        % partial of the measure wrt the landmark position:
        h_pi=H_Pi(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,S0(vis(i),1),S0(vis(i),2),S0(vis(i),3),s1,s2,s3,sc(1),sc(2),sc(3),x0,y0,z0);
        
        keep_vis_index_covmat=[ns*match(i,1)-(ns-1)+[0:(ns-1)]]';
          
        Hi=[h_xi,h_pi];
        Pi=[Prr,Prm(:,keep_vis_index_covmat);Prm(:,keep_vis_index_covmat)',Pmm(keep_vis_index_covmat,keep_vis_index_covmat)];
        Zi=Hi*Pi*Hi';
        meas_uncert(i)=det(Zi);
        end
        
        [~,order]=sort(meas_uncert);
        
        % max 10 updates each time:
        to_correct=vis(order(1:min(10,end)));

end