%% initialize variables
initialize_params;



% 
% pix_coord=double(measures{loop}); % convert to double precision
% y=[pix_coord([1,2],:);pix_coord(3,:)-pix_coord(1,:)];
% yn.m=size(y,2);
% yn.z=reshape(y,[],1);
% yn.feats=feats{loop};


% load('meas.mat');
% 
% %% preprocessing measures:
% for i=1:size(measures,1)
%     measures{i}([1,3],:)=measures{i}([1,3],:)-cam_params.hpix/2;
%     measures{i}([2,4],:)=-measures{i}([2,4],:)+cam_params.vpix/2;
% end

% %%
% tf=Nmax*dt;
% 
%     % setting tolerance for numerical integration:
%         options=odeset('AbsTol',1e-10,'RelTol',1e-12);
%     % ode integration:
%         [T,Y]=ode113(@(t,y) process(t,y,params),[0,tf],x0,options);
%         
%         for i=1:size(Y,1)
%             [y]=sim_measure(xn,fv,fv2,cam_params,MASK);
% %             
% %             xr=Y(i,1);
% %             yr=Y(i,2);
% %             zr=Y(i,3);
% %             theta=Y(i,9);
% %             s=Y(i,14:16)';
% %             s1=s(1);
% %             s2=s(2);
% %             s3=s(3);
% % 
% %         % define skew symmetric matrix as cross(s,.).
% %             skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
% %         % define DCM that brings from chaser to target frame D:
% %             D=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;
% %         % D' is DCM from target to chaser:
% %             D=D';
% %             C_BI=quat2dcm(qc);
% %             C_LI=[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
% %             C_BL=C_BI*C_LI';
% %         % S is the array containing all points in the complete model:
% %             P_i=C_BL'*D*fv2.Points'+[xr;yr;zr];
% %             fv_new2=triangulation(fv2.ConnectivityList,P_i');
% %             trisurf(fv_new2)
% %             hold on
%         end
% 


%%
% cam=cam_params;
%     s=xn(14:16);
%     s1=s(1);
%     s2=s(2);
%     s3=s(3);
%     r=xn(7);
%     f=xn(9);
%     
%     
%     
%     sc=xn(22:end);
%             sc1=sc(1);
%             sc2=sc(2);
%             sc3=sc(3);
%             
%             
%             
%             skew_sc=[0,-sc3,sc2;sc3,0,-sc1;-sc2,sc1,0];
%         % define DCM that brings from I to C frame:
%             C_CI=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;
%             C_LI=[cos(f),sin(f),0;-sin(f),cos(f),0;0,0,1];
%             C_CL=C_CI*C_LI';
%             
%     
%     
%     skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
% % define DCM that brings from chaser to target frame D:
%     C_TC=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;
% 
%     P_i=(C_TC'*fv.Points'+C_CL*xn(1:3)); % (Pi-T)+(T-C)
% 
%     % definition of the rotated triangulations:
%     fv_new=triangulation(fv.ConnectivityList,P_i');
% 
%     xi=P_i(1,:);
%     yi=P_i(2,:);
%     zi=P_i(3,:);
% 
%     % simula l'acquisizione in camera (pinhole camera model)+ disparity(stereo)
%     h=[cam.u0-cam.hpix*cam.foc/cam.Hf*xi./zi;cam.v0-cam.vpix*cam.foc/cam.Vf*yi./zi;cam.foc*cam.b./zi];
%     % serve perchè in stereovision la disparità maggiore è data ai punti più
%     % vicini. (solo per visualizzare i plot).
%     transformedh=[h(1:2,:);1./h(3,:)];
% 
%     % definition of the rotated triangulations:
%     fv_cam=triangulation(fv.ConnectivityList,transformedh');
% 
%     trisurf(fv_cam);
% 
%     figure(2)
%     title('Simulated Image');
%     trisurf(fv_cam);
% 
%     hold on
%     plot3(0,0,0,'ro','LineWidth',2);
%     axis equal
%     campos([0,0,0]);
%     % view([0,90]);
%     xlim([-cam.hpix/2,cam.hpix/2]);
%     ylim([-cam.vpix/2,cam.vpix/2]);
%     ax = gca;
%     % ax.YDir = 'reverse';
%     % ax.XDir = 'reverse';
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');


%%


%% Start Loop:
x_hist=[];
r_hist=[];
lmkinfo.counter_meas=[];
lmkinfo.counter_prop=[];
feats_list=single(zeros(0,64));

for loop=1:Nmax-start_frame+1
    loop
    %% Simulator step:
    % obtain new state and new estimated measure after dt:
    % x0, xn, yn are 'reality':
%         [xn,yn]=simulate_next_step(x0,dt,fv,fv2,cam_params,qc,MASK);
        xn=sim_states(x0,dt,params);

%% image processing part:


        [yn]=process_images(renderer,loop+start_frame-1,proc_method);



    %% Filter propagator step:
    % X0, Xn, Yn are estimates and S0 is set of available landmarks, Sn is 
    % the set of points which are estimated to be visible in FOV:
        [Xn,Yn,Prrn,Prmn,Pmmn,lmkinfo]=propagate_next_step(X0,Prr0,Prm0,Pmm0,Q,dt,S0,cam_params,R,lmkinfo,params);        
        


%         meas_act=reshape(yn.z,3,[])';
%         scatter3(meas_act(:,1),meas_act(:,2),meas_act(:,3));
%         meas_est=reshape(Yn.h,3,[])';
%         hold on
%         scatter3(meas_est(:,1),meas_est(:,2),meas_est(:,3));


        compatibility=compute_compatibility (Yn, yn);
        H= JCBB_test(Yn,yn,compatibility);
%         sum(find(H))
%         [~,i, j] = find(H)

        [Prrn1,Prmn1,Pmmn1,Xn1,Sn1,lmkinfo,feats_list]=update_step(Xn,Prrn,Prmn,Pmmn,yn,Yn,S0,H,Opt,lmkinfo,feats_list);
   
%%

        new = find((H == 0) & (compatibility.AL == 0));
        new=new(1:min(Opt.Nnew(min(loop,2)),length(new)));
        
%         close all
%         figure()

%         meas_est=reshape(Yn.h,3,[])';
%         hold on
%         scatter(meas_est(:,1),meas_est(:,2));
%         
        
%         plot(meas_act(:,3))
%         hold on
%         plot(meas_est(:,3))
        
        
        if size(S0,1)<max_landmarks
            [S_new,Prm_new,Pmm_new,lmkinfo,feats_list]=add_points(Xn1,yn,Prrn1,Prmn1,Pmmn1,R,cam_params,new,lmkinfo,feats_list);
        else
            S_new=[];
            Prm_new=Prmn1;
            Pmm_new=Pmmn1;
        end
        


%%







        
        
        [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn1,Prrn1,Prm_new,Pmm_new,S_new,Sn1);
        [S0,Prm0,Pmm0,lmkinfo,feats_list]=remove_points(S0,Prm0,Pmm0,lmkinfo,feats_list);
        
        
        
               [simulator_history,state_history,landmarks_map,covariance_mat,visibility_check]=Store_Iteration(xn,yn,...
        Xn,S0,Yn,Prrn,Prmn,Pmmn,...
        Prrn1,Prmn1,Pmmn1,Xn1,Sn1,...
        Prm_new,Pmm_new,S_new,...
        simulator_history,state_history,landmarks_map,covariance_mat,loop,visibility_check);

        
        %% live plot
        
        
        live_processing;
        
        
        

end

%% simple kalman filter application:

% 
%     [xn,yn]=simulate_next_step(x0,dt,fv,fv2,cam_params,qc);
% %% Filter propagator step:
% % X0, Xn, Yn are estimates and S0 is set of available landmarks, Sn is 
% % the set of points which are estimated to be visible in FOV:
% %     X0=x0;
% %     [Xn,Yn,Prrn,Prmn,Pmmn]=propagate_next_step(X0,Prr0,Prm0,Pmm0,Q,dt,S0,cam_params,qc,R);        
% %     compatibility=compute_compatibility (Yn, yn);
% %     H= JCBB_test(Yn,yn,compatibility);
% % 
% %     [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=update_step(Xn,Prrn,Prmn,Pmmn,yn,Yn,S0,H,Opt);
% % 
% %     new = find((H == 0) & (compatibility.AL == 0));
% %     new=new(1:10);
% %     [S_new,Prm_new,Pmm_new]=add_points(Xn1,yn,Prrn1,Prmn1,Pmmn1,R,cam_params,qc,new); 
% %     [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn1,Prrn1,Prm_new,Pmm_new,S_new,Sn1);        
% spr=0.1;    
% 
%     Pmm0=1*eye(15);
%     Prm0=zeros(6,15);
% 
%     close all
% % scatter3(S0(:,1),S0(:,2),S0(:,3))
% % hold on
% % scatter3(fv2.Points(1:5,1),fv2.Points(1:5,2),fv2.Points(1:5,3))
% % 
% 
% 
% for loop=1:Nmax
%     loop
% 
%     %% Simulator step:
%     % obtain new state and new estimated measure after dt:
%     % x0, xn, yn are 'reality':
%         [xn,yn]=simulate_next_step(x0,dt,fv,fv2,cam_params,qc);
%     %% Filter propagator step:
%     % X0, Xn, Yn are estimates and S0 is set of available landmarks, Sn is 
%     % the set of points which are estimated to be visible in FOV:
%         [Xn,Yn,Prrn,Prmn,Pmmn]=propagate_next_step(X0,Prr0,Prm0,Pmm0,Q,dt,S0,cam_params,qc,R);
%         
%         [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=update_step_noslam(Xn,Prrn,Prmn,Pmmn,yn,Yn,S0);
%         
%         S_new=[];
%         Prm_new=Prmn1;
%         Pmm_new=Pmmn1;
%         
% %         norm(xn(1:3)-Xn1(1:3))
% %         norm(xn(4:6)-Xn1(4:6))
% %         norm(xn(7:end)-Xn1(7:end))
%         x_hist=[x_hist;Xn1(1:6)];
%         r_hist=[r_hist;xn(1:6)];
%         
% %         plot(meas_act(:,3))
% %         hold on
% %         plot(meas_est(:,3))
%         
%         [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn1,Prrn1,Prm_new,Pmm_new,S_new,Sn1);        
% 
%         scatter3(S0(:,1),S0(:,2),S0(:,3))
%         hold on;
% end
% 
% state_ekf=reshape(x_hist,6,[]);
% state_real=reshape(r_hist,6,[]);
% 
% figure()
% plot(vecnorm(state_ekf(1,:)-state_real(1,:),2,1));
% hold on
% plot(vecnorm(state_ekf(2,:)-state_real(2,:),2,1));
% plot(vecnorm(state_ekf(3,:)-state_real(3,:),2,1));
% 
% figure()
% plot(vecnorm(state_ekf(4,:)-state_real(4,:),2,1));
% hold on
% plot(vecnorm(state_ekf(5,:)-state_real(5,:),2,1));
% plot(vecnorm(state_ekf(6,:)-state_real(6,:),2,1));
%  


%                
%         
%         
%         for i=1:size(Yn,2)
%             dz=Yn(:,1)-yn;
%             
%             
%             
%         end
%         
%         
%         
%         
%         
%         
%         
%         
%         
% %         scatter(yn(1,:),yn(2,:))
% %         hold on
% %         scatter(Yn(1,:),Yn(2,:))
% 
%     %% Matching of measures data aggregation algorithm:
%     % matching between estimated visible points and real measures:
%         tic;
% %       [match,Hn,Zn]=match_measures(Xn,yn,Yn,Prrn,Prmn,Pmmn,R,S0,vis,b,foc,qc);
%         [match,Hn,Zn]=automatch(index_meas,vis,map_visibility,Xn,S0,Prrn,Prmn,Pmmn,R,b,foc,qc);
%         t_match=toc;
% 
%     %% Filter update:
%     % update states, landmarks and covariances using matched measures:
%         tic;
%         [to_correct]=SelectLandmarksToCorrect(Xn,Prrn,Prmn,Pmmn,vis,qc,cam_params);
%         [matched]=MatchingRANSAC();
% %         [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=update_step(Xn,Prrn,Prmn,Pmmn,Zn,Hn,yn,Yn,S0,vis,match);
%         [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=correction_one_by_one(b,foc,Xn,Prrn,Prmn,Pmmn,yn,Yn,S0,R,vis,match,qc);
%         t_update=toc;
% %         eig(full([Prrn1,Prmn1;Prmn1',Pmmn1]))
% 
%     %% New landmarks initialization:
%     % initialize new landmarks and their covariance:
%         tic;
%         % sottraggo le misure matchate sia dall'index_meas che dalla misura
%         index_meas(match(:,2))=[];
%         yn(:,match(:,2))=[];
%         non_new_feats=find(ismember(index_meas,map_visibility));
%         index_meas(non_new_feats)=[];
%         yn(:,non_new_feats)=[];
%         n=min(1,size(yn,2));
%         map_visibility=[map_visibility;index_meas(1:n)];
%         
% %         [S_new,Prm_new,Pmm_new]=add_points(Xn1,yn,b,foc,Prrn1,Prmn1,Pmmn1,R,qc);
%         [S_new,Prm_new,Pmm_new]=initialize_landmarks_one_by_one(Xn1,yn,b,foc,Prrn1,Prmn1,Pmmn1,R,qc);
%         t_init=toc;
% %         eig(full([Prrn,Prm_new;Prm_new',Pmm_new]))
%         
%         
%     %% Iteration storage:
%         [t_storage,times_history,simulator_history,state_history,landmarks_map,covariance_mat,visibility_check]=Store_Iteration(t_sim,t_match,t_update,t_init,xn,yn,...
%         Xn,S0,Yn,Prrn,Prmn,Pmmn,vis,...
%         match,Hn,Zn,...
%         Prrn1,Prmn1,Pmmn1,Xn1,Sn1,...
%         Prm_new,Pmm_new,S_new,...
%         times_history,simulator_history,state_history,landmarks_map,covariance_mat,visibility_check,loop);
% 
%     %% reinitialization:
%     % reset initial values for next iteration:
%         [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn1,Prrn1,Prm_new,Pmm_new,S_new,Sn1);
%     
%     %% message information:
%         fprintf('\t Iteration nr. %d  \n', loop);
%         fprintf('------------------------------------------\n');
%         fprintf('->\t Simulation step          : %.3f  [s]\n', t_sim);
%         fprintf('->\t Filter propagation step  : %.3f  [s]\n', t_prop);
%         fprintf('->\t Matching measure step    : %.3f  [s]\n', t_match);
%         fprintf('->\t Filter update step       : %.3f  [s]\n', t_update);
%         fprintf('->\t Landmarks Initialization : %.3f  [s]\n', t_init);
%         fprintf('->\t Iteration storage        : %.3f  [s]\n', t_storage);
%         fprintf('------------------------------------------\n');
%         fprintf('\t Tot.                     : %.3f  [s]\n', t_storage+t_sim+t_prop+t_match+t_update+t_init);
% end

%% Postprocessing

% Nmax=loop;
% validation_script;
% plotting_script;

