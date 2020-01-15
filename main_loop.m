

%% initialize variables
initialize_params;

%% Start Loop:
x_hist=[];
r_hist=[];
lmkinfo.counter_meas=[];
lmkinfo.counter_prop=[];
feats_list=single(zeros(0,64));

loop=1;
[y0,frameLeftGray0]=process_images(renderer,loop+start_frame-2,proc_method,cam_params);

for loop=1:147
    
    xn=sim_states(x0,dt,params);
    [yn,frameLeftGrayn]=process_images(renderer,loop+start_frame-1,proc_method,cam_params);
    indexPairs = matchFeatures(y0.feats,yn.feats);
    
    [Xn]=prop_states(X0,dt,params);
    [hn]=new_meas(y0,X0,Xn,cam_params,indexPairs);
    %     h0=p0(indexPairs(:,1));
    
    mn=reshape(hn',[],1);
    %
    p0=reshape(y0.z,3,[])';
    pn=reshape(yn.z,3,[])';
    figure
    showMatchedFeatures(frameLeftGray0,frameLeftGrayn,p0(indexPairs(:,1),1:2),pn(indexPairs(:,2),1:2));
    
    
    % covariance computation step:
    Phi0=STM(X0,dt);
    P1=Phi0*Prr0*Phi0'+Q;
    [A,B,C]=matrix_construction(X0,Xn,p0(indexPairs(:,1),:),cam_params);
    
    RkCell = repmat({R}, 1, size(indexPairs,1));
    Rk = blkdiag(RkCell{:});
    
    Zn= A*P1*A' + B*P0*B' + A*Phi0*P0*B' + B*P0*Phi0*A' + C*Rk*C' +Rk;
    
    
    
    
    %     % debug matching across time instants:
    
    %     p0=reshape(y0.z,3,[])';
    %     pn=reshape(yn.z,3,[])';
    
    
    %%
    
    
    
    [S0,Prm0,Pmm0]=add_points(X0,y0,Prr0,Prm0,Pmm0,R,cam_params,indexPairs);
    
    [Xn]=prop_states(X0,dt,params);
    
    % propagate covariance:
    % covariance propagation step:
    
    % state only covariance propagation:
    Prrn=Phi0*Prr0*Phi0'+Q;
    % mixed state landmarks covariance propagation:
    Prmn=Phi0*Prm0;
    % the covariance propagation is independent of the collected points and
    % does not influence their pure covariance Pmm, just the mixed state terms.
    Pmmn=Pmm0;
    
    % estimate new measures at estimated new state:
    % estimate of measures and visibility (inside FOV):
    [Yn] = prop_measures(Xn,S0,Prrn,Prmn,Pmmn,cam_params,R);
    
    
    [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=update_step(Xn,Prrn,Prmn,Pmmn,yn,Yn,S0,indexPairs);
    
    % estimated state:
    y0=yn;
    X0=Xn1;
    x0=xn;
    
    S0=zeros(0,3);
    Prr0=0.5*(Prrn1+Prrn1');
    Prm0=sparse(zeros(14,3*size(S0,1)));
    Pmm0=sparse(zeros(3*size(S0,1)));
    
    
    
    live_processing;
    
    
end
%%

%
% for loop=1:Nmax-start_frame+1
%     loop
%     %% Simulator step:
%     % obtain new state and new estimated measure after dt:
%     % x0, xn, yn are 'reality':
% %         [xn,yn]=simulate_next_step(x0,dt,fv,fv2,cam_params,qc,MASK);
%         xn=sim_states(x0,dt,params);
%
% %% image processing part:
%
%
%         [yn]=process_images(renderer,loop+start_frame-1,proc_method,cam_params);
%
%
%
%
%
%
%
%
%
%
%
%
%         %% matching between old and new measures:
%
%         [mid]=match_measures(y0,yn);
%
%
%         % initialize matched states
%         if size(S0,1)<max_landmarks
%             [S_new,Prm_new,Pmm_new,lmkinfo,feats_list]=add_points(X0,yn,Prr0,Prm0,Pmm0,R,cam_params,mid,lmkinfo,feats_list);
%         else
%             S_new=[];
%             Prm_new=Prmn1;
%             Pmm_new=Pmmn1;
%         end
%
%     % Filter propagator step:
%     % X0, Xn, Yn are estimates and S0 is set of available landmarks, Sn is
%     % the set of points which are estimated to be visible in FOV:
%
%     % step to propagate the states inside the filter:
%         [Xn]=prop_states(X0,dt,params);
%
%     % propagate covariance:
%     % covariance propagation step:
%         Phi0=STM(X0,dt);
%     % state only covariance propagation:
%         Prrn=Phi0*Prr0*Phi0'+Q0;
%     % mixed state landmarks covariance propagation:
%         Prmn=Phi0*Prm0;
%     % the covariance propagation is independent of the collected points and
%     % does not influence their pure covariance Pmm, just the mixed state terms.
%         Pmmn=Pmm0;
%
%     % estimate new measures at estimated new state:
%     % estimate of measures and visibility (inside FOV):
%         [Yn,lmkinfo]=prop_measures(Xn,S0,Prrn,Prmn,Pmmn,cam_params,R,lmkinfo);
%
%
%
%
%
%
%
%
%
%
%         [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn1,Prrn1,Prm_new,Pmm_new,S_new,Sn1);
%         [S0,Prm0,Pmm0,lmkinfo,feats_list]=remove_points(S0,Prm0,Pmm0,lmkinfo,feats_list);
%
%         [simulator_history,state_history,landmarks_map,covariance_mat,visibility_check]=Store_Iteration(xn,yn,...
%         Xn,S0,Yn,Prrn,Prmn,Pmmn,...
%         Prrn1,Prmn1,Pmmn1,Xn1,Sn1,...
%         Prm_new,Pmm_new,S_new,...
%         simulator_history,state_history,landmarks_map,covariance_mat,loop,visibility_check);
%
%
%         %% live plot
%
%         live_processing;
%
%
% end
%
% %% simple kalman filter application:
%
% %
% %     [xn,yn]=simulate_next_step(x0,dt,fv,fv2,cam_params,qc);
% % %% Filter propagator step:
% % % X0, Xn, Yn are estimates and S0 is set of available landmarks, Sn is
% % % the set of points which are estimated to be visible in FOV:
% % %     X0=x0;
% % %     [Xn,Yn,Prrn,Prmn,Pmmn]=propagate_next_step(X0,Prr0,Prm0,Pmm0,Q,dt,S0,cam_params,qc,R);
% % %     compatibility=compute_compatibility (Yn, yn);
% % %     H= JCBB_test(Yn,yn,compatibility);
% % %
% % %     [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=update_step(Xn,Prrn,Prmn,Pmmn,yn,Yn,S0,H,Opt);
% % %
% % %     new = find((H == 0) & (compatibility.AL == 0));
% % %     new=new(1:10);
% % %     [S_new,Prm_new,Pmm_new]=add_points(Xn1,yn,Prrn1,Prmn1,Pmmn1,R,cam_params,qc,new);
% % %     [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn1,Prrn1,Prm_new,Pmm_new,S_new,Sn1);
% % spr=0.1;
% %
% %     Pmm0=1*eye(15);
% %     Prm0=zeros(6,15);
% %
% %     close all
% % % scatter3(S0(:,1),S0(:,2),S0(:,3))
% % % hold on
% % % scatter3(fv2.Points(1:5,1),fv2.Points(1:5,2),fv2.Points(1:5,3))
% % %
% %
% %
% % for loop=1:Nmax
% %     loop
% %
% %     %% Simulator step:
% %     % obtain new state and new estimated measure after dt:
% %     % x0, xn, yn are 'reality':
% %         [xn,yn]=simulate_next_step(x0,dt,fv,fv2,cam_params,qc);
% %     %% Filter propagator step:
% %     % X0, Xn, Yn are estimates and S0 is set of available landmarks, Sn is
% %     % the set of points which are estimated to be visible in FOV:
% %         [Xn,Yn,Prrn,Prmn,Pmmn]=propagate_next_step(X0,Prr0,Prm0,Pmm0,Q,dt,S0,cam_params,qc,R);
% %
% %         [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=update_step_noslam(Xn,Prrn,Prmn,Pmmn,yn,Yn,S0);
% %
% %         S_new=[];
% %         Prm_new=Prmn1;
% %         Pmm_new=Pmmn1;
% %
% % %         norm(xn(1:3)-Xn1(1:3))
% % %         norm(xn(4:6)-Xn1(4:6))
% % %         norm(xn(7:end)-Xn1(7:end))
% %         x_hist=[x_hist;Xn1(1:6)];
% %         r_hist=[r_hist;xn(1:6)];
% %
% % %         plot(meas_act(:,3))
% % %         hold on
% % %         plot(meas_est(:,3))
% %
% %         [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn1,Prrn1,Prm_new,Pmm_new,S_new,Sn1);
% %
% %         scatter3(S0(:,1),S0(:,2),S0(:,3))
% %         hold on;
% % end
% %
% % state_ekf=reshape(x_hist,6,[]);
% % state_real=reshape(r_hist,6,[]);
% %
% % figure()
% % plot(vecnorm(state_ekf(1,:)-state_real(1,:),2,1));
% % hold on
% % plot(vecnorm(state_ekf(2,:)-state_real(2,:),2,1));
% % plot(vecnorm(state_ekf(3,:)-state_real(3,:),2,1));
% %
% % figure()
% % plot(vecnorm(state_ekf(4,:)-state_real(4,:),2,1));
% % hold on
% % plot(vecnorm(state_ekf(5,:)-state_real(5,:),2,1));
% % plot(vecnorm(state_ekf(6,:)-state_real(6,:),2,1));
% %
%
%
% %
% %
% %
% %         for i=1:size(Yn,2)
% %             dz=Yn(:,1)-yn;
% %
% %
% %
% %         end
% %
% %
% %
% %
% %
% %
% %
% %
% %
% % %         scatter(yn(1,:),yn(2,:))
% % %         hold on
% % %         scatter(Yn(1,:),Yn(2,:))
% %
% %     %% Matching of measures data aggregation algorithm:
% %     % matching between estimated visible points and real measures:
% %         tic;
% % %       [match,Hn,Zn]=match_measures(Xn,yn,Yn,Prrn,Prmn,Pmmn,R,S0,vis,b,foc,qc);
% %         [match,Hn,Zn]=automatch(index_meas,vis,map_visibility,Xn,S0,Prrn,Prmn,Pmmn,R,b,foc,qc);
% %         t_match=toc;
% %
% %     %% Filter update:
% %     % update states, landmarks and covariances using matched measures:
% %         tic;
% %         [to_correct]=SelectLandmarksToCorrect(Xn,Prrn,Prmn,Pmmn,vis,qc,cam_params);
% %         [matched]=MatchingRANSAC();
% % %         [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=update_step(Xn,Prrn,Prmn,Pmmn,Zn,Hn,yn,Yn,S0,vis,match);
% %         [Prrn1,Prmn1,Pmmn1,Xn1,Sn1]=correction_one_by_one(b,foc,Xn,Prrn,Prmn,Pmmn,yn,Yn,S0,R,vis,match,qc);
% %         t_update=toc;
% % %         eig(full([Prrn1,Prmn1;Prmn1',Pmmn1]))
% %
% %     %% New landmarks initialization:
% %     % initialize new landmarks and their covariance:
% %         tic;
% %         % sottraggo le misure matchate sia dall'index_meas che dalla misura
% %         index_meas(match(:,2))=[];
% %         yn(:,match(:,2))=[];
% %         non_new_feats=find(ismember(index_meas,map_visibility));
% %         index_meas(non_new_feats)=[];
% %         yn(:,non_new_feats)=[];
% %         n=min(1,size(yn,2));
% %         map_visibility=[map_visibility;index_meas(1:n)];
% %
% % %         [S_new,Prm_new,Pmm_new]=add_points(Xn1,yn,b,foc,Prrn1,Prmn1,Pmmn1,R,qc);
% %         [S_new,Prm_new,Pmm_new]=initialize_landmarks_one_by_one(Xn1,yn,b,foc,Prrn1,Prmn1,Pmmn1,R,qc);
% %         t_init=toc;
% % %         eig(full([Prrn,Prm_new;Prm_new',Pmm_new]))
% %
% %
% %     %% Iteration storage:
% %         [t_storage,times_history,simulator_history,state_history,landmarks_map,covariance_mat,visibility_check]=Store_Iteration(t_sim,t_match,t_update,t_init,xn,yn,...
% %         Xn,S0,Yn,Prrn,Prmn,Pmmn,vis,...
% %         match,Hn,Zn,...
% %         Prrn1,Prmn1,Pmmn1,Xn1,Sn1,...
% %         Prm_new,Pmm_new,S_new,...
% %         times_history,simulator_history,state_history,landmarks_map,covariance_mat,visibility_check,loop);
% %
% %     %% reinitialization:
% %     % reset initial values for next iteration:
% %         [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn1,Prrn1,Prm_new,Pmm_new,S_new,Sn1);
% %
% %     %% message information:
% %         fprintf('\t Iteration nr. %d  \n', loop);
% %         fprintf('------------------------------------------\n');
% %         fprintf('->\t Simulation step          : %.3f  [s]\n', t_sim);
% %         fprintf('->\t Filter propagation step  : %.3f  [s]\n', t_prop);
% %         fprintf('->\t Matching measure step    : %.3f  [s]\n', t_match);
% %         fprintf('->\t Filter update step       : %.3f  [s]\n', t_update);
% %         fprintf('->\t Landmarks Initialization : %.3f  [s]\n', t_init);
% %         fprintf('->\t Iteration storage        : %.3f  [s]\n', t_storage);
% %         fprintf('------------------------------------------\n');
% %         fprintf('\t Tot.                     : %.3f  [s]\n', t_storage+t_sim+t_prop+t_match+t_update+t_init);
% % end
%
% %% Postprocessing
%
% % Nmax=loop;
% % validation_script;
% % plotting_script;

