function [y]=sim_measure(X,fv,fv2,cam_params,MASK)
    % This function simulates a time step in the simulator.
    % Inputs:
    % X    :  state of simulator after dt;
    % fv   :  triangulation object of detailed model;
    % fv2  :  triangulation object of reduced model;
    % npix :  number of pixels on vertical side of sensor;
    % Hf   :  horizontal size of sensor [m];
    % Vf   :  vertical size of sensor [m];
    % foc  :  focal length for camera [m];
    % b    :  baseline between cameras along y [m];
    %
    % Outputs:
    % y    : vector of size [nh,m] of m measurements of size nh each;

    %% conversion of points from target to chaser frame:
    % fv2 for points, fv for triangles
    % retrieve position and relative orientation from current simulator state
    % vector:
        x0=X(1);
        y0=X(2);
        z0=X(3);
        theta=X(9);
        s=X(14:16);
        s1=s(1);
        s2=s(2);
        s3=s(3);

    % define skew symmetric matrix as cross(s,.).
        skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
    % define DCM that brings from chaser to target frame D:
        D=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;
    % D' is DCM from target to chaser:
        D=D';
    % S is the array containing all points in the complete model:
        S=fv.Points;

    % definition of the DCM that rotates from LVLH to chaser body frame:
%     sc=[sc1;sc2;sc3];
%         skew_sc=[0,-sc(3),sc(2);sc(3),0,-sc(1);-sc(2),sc(1),0];
%         Dc=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;
%         C_BL=simplify(Dc);
%         C_BI=quat2dcm(qc);
%         C_LI=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
%         C_BL=C_BI*C_LI';
%         qc=dcm2quat(C_BL);
%         sc=qc(2:end)/(1+qc(1));
        sc=X(22:24);
        skew_sc=[0,-sc(3),sc(2);sc(3),0,-sc(1);-sc(2),sc(1),0];
        Dc=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;
        

        C_BI=(Dc);
        C_LI=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
        C_BL=C_BI*C_LI';

    % definition of points of the complete and reduced model in the camera
    % frame (i.e. aligned with x axis of chaser):
        P_i=D*S'+C_BL*[x0;y0;z0];
        P_i2=D*fv2.Points'+C_BL*[x0;y0;z0];

    % definition of the rotated triangulations:
        fv_new=triangulation(fv.ConnectivityList,P_i');
        fv_new2=triangulation(fv2.ConnectivityList,P_i2');

    % extraction of the first, second, and third point of each triangle from
    % the detailed model. extracted as arrays for GPU implementation:
        P0=fv_new2.Points(fv_new2.ConnectivityList(:,1),:);
        P1=fv_new2.Points(fv_new2.ConnectivityList(:,2),:);
        P2=fv_new2.Points(fv_new2.ConnectivityList(:,3),:);

    % definition of matrix of directions from 0 (chaser center of mass) to each
    % of the reduced model points. Also definition of directions from second
    % camera (placed at -b along y axis):
        mat_D=fv_new2.Points;
        mat_D2=mat_D-[cam_params.b,0,0];
        pnorm=vecnorm(mat_D,2,2);
        pnorm2=vecnorm(mat_D2,2,2);
        mat_D=mat_D./pnorm;
        mat_D2=mat_D2./pnorm2;

    % conversion in gpuarrays for gpu implementation:
        gD=gpuArray(mat_D);
        gD2=gpuArray(mat_D2);

    %% Ray tracing step:
    % gpu array function that detects intersection of infinite rays with
    % defined triangles for bot left and right camera:
        [dist, flag] = arrayfun(@rayTriGPU, P0(:,1)', P0(:,2)', P0(:,3)', ...
        P1(:,1)', P1(:,2)', P1(:,3)', ...
        P2(:,1)', P2(:,2)', P2(:,3)', ...
        0, 0, 0, ...
        gD(:,1),gD(:,2),gD(:,3));
        [dist2, flag2] = arrayfun(@rayTriGPU, P0(:,1)', P0(:,2)', P0(:,3)', ...
        P1(:,1)', P1(:,2)', P1(:,3)', ...
        P2(:,1)', P2(:,2)', P2(:,3)', ...
        cam_params.b, 0, 0, ...
        gD2(:,1),gD2(:,2),gD2(:,3));

    %% processing of intersections for both cameras:
    % conversion to normal array (from gpu):
        distances = gather(dist);
    % substitute NaN with inf in array of distances from origin of rays:
        distances(~isfinite(distances)) = inf;
        distances(find(MASK)) = inf;
    % for each ray extract the minimum distance of intersection:
        minint=min(distances,[],2);

    % conversion to normal array (from gpu):
        distances2 = gather(dist2);
    % substitute NaN with inf in array of distances from origin of rays. Nan
    % values appear if no intersection is found (e.g. point is perfectly
    % visible and is not in front or behind anything).
        distances2(~isfinite(distances2)) = inf;
        distances2(find(MASK)) = inf;
    % for each ray extract the minimum distance of intersection:
        minint2=min(distances2,[],2);

    %% visibility check:
    % in order for the feature to be visible it must be visible for both
    % cameras. Means that the minimum distance from intersection must be larger
    % than the distance from the camera (i.e. feature is visible but in front
    % of one of the triangular surfaces:
        vis_index=find(pnorm2<minint2 & pnorm<minint);

        points=fv_new2.Points(pnorm2<minint2 & pnorm<minint,:);
        
        %%%%% WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         points=fv_new2.Points;
%         vis_index=linspace(1,size(points,1),size(points,1))';

    %% measurements conversion to pixels:
        P_i=points';
        xi=P_i(1,:);
        yi=P_i(2,:);
        zi=P_i(3,:);
    % stereo conversion measurement equation:
        h=[cam_params.u0+cam_params.alpha_u*xi./zi;cam_params.v0+cam_params.alpha_v*yi./zi;cam_params.alpha_u*cam_params.b./zi];
        scatter(h(1,:),h(2,:));
        xlim([-cam_params.hpix/2,cam_params.hpix/2]);
        ylim([-cam_params.vpix/2,cam_params.vpix/2]);

%         H_fun(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,fv2.Points(1,1),fv2.Points(1,2),fv2.Points(1,3),s1,s2,s3,sc(1),sc(2),sc(3),cam_params.u0,cam_params.v0,x0,y0,z0)
        %         scatter(h(1,:),h(2,:))
    % points are registered by pixels if they are in FOV of camera (between +
    % and - the size of the sensor):
        vis=abs(h(1,:))<cam_params.hpix/2 & abs(h(2,:))<cam_params.vpix/2;
        in_FOV_index=vis_index(vis);
        vis_meas=h(:,vis>0);
%         vis_meas(1:2,:)=fix(vis_meas(1:2,:));
        
%         y=[vis_meas];
        
        
        y.m=size(vis_meas,2);
        y.z=reshape(vis_meas,[],1);
%         scatter(vis_meas(1,1:10),vis_meas(2,1:10))
        y.vis=in_FOV_index;       

end