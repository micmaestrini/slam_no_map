

%% initialize variables
initialize_params;

%% Start Loop:
x_hist=[];
r_hist=[];
lmkinfo.counter_meas=[];
lmkinfo.counter_prop=[];
feats_list=single(zeros(0,64));

loop=1;
[frameLeftGray0,filt_match10,filt_match20,feats_HK0]=process_images(renderer,loop+start_frame-2,proc_method,cam_params);

for loop=1:147
    
    xn=sim_states(x0,dt,params);
    
    
    [frameLeftGrayn,filt_match1n,filt_match2n, feats_HKn]=process_images(renderer,loop+start_frame-1,proc_method,cam_params);
    
    
    indexPairs = matchFeatures(feats_HK0,feats_HKn,'Unique',true);
    
    % match between time frames including outliers
    outliern1=filt_match1n(indexPairs(:,2),:);
    outliern2=filt_match2n(indexPairs(:,2),:);
    outlier01=filt_match10(indexPairs(:,1),:);
    outlier02=filt_match20(indexPairs(:,1),:);
    
    
    [tform, inliern1, inlier01] = estimateGeometricTransform(...
        outliern1, outlier01, 'similarity','MaxDistance',5);
    
    
    [~,d]=ismember(inliern1.Location,outliern1.Location);
    
    xq0=outlier01.Location(d(:,1));
    yq0=outlier01.Location(d(:,2));
    xqr0=outlier02.Location(d(:,1));
    vq0=xqr0-xq0;
    measures0=[xq0-cam_params.hpix/2,-yq0+cam_params.vpix/2,vq0];
    y0.m=size(measures0,1);
    y0.z=double(reshape(measures0',[],1));
    
    
    
    xqn=outliern1.Location(d(:,1));
    yqn=outliern1.Location(d(:,2));
    xqrn=outliern2.Location(d(:,1));
    vqn=xqrn-xqn;
    measuresn=[xqn-cam_params.hpix/2,-yqn+cam_params.vpix/2,vqn];
    yn.m=size(measuresn,1);
    yn.z=double(reshape(measuresn',[],1)); % cast to double to use in functions
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    showMatchedFeatures(frameLeftGray0,frameLeftGrayn,outlier01,outliern1);
    figure;
    showMatchedFeatures(frameLeftGray0,frameLeftGrayn,inlier01,inliern1);
    title('Matching points (inliers only)');
    legend('ptsOriginal','ptsDistorted');
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Xn]=prop_states(X0,dt,params);
    %
    p0=reshape(y0.z,3,[])';
    
    %     mea_0=reshape(p0(indexPairs(:,1),:)',[],1);
    %     mea_n=reshape(pn(indexPairs(:,2),:)',[],1);
    
    % covariance computation step:
    Phi0=STM(X0,dt);
    Prrn=Phi0*Prr0*Phi0'+Q;
    [A,B,C,h_x]=matrix_construction(X0,Xn,p0,cam_params);
    
    RkCell = repmat({R}, 1, yn.m);
    Rk = blkdiag(RkCell{:});
    
    Zn= A*Prrn*A' + B*Prr0*B' + A*Phi0*Prr0*B' + B*Prr0*Phi0'*A' + C*Rk*C' +Rk;
    Zn= A*Prrn*A'+ B*Prr0*B'+ C*Rk*C'+ Rk;
    zn=yn.z-h_x;
    
    % update step
    staten=Xn([1:6,11:18]);
    params1=Xn(7:10);
    params2=Xn(19:24);
    
    K=Prrn*A'*inv(Zn);
    state1=staten+K*zn;
    Prr1=Prrn-K*A*Prrn;
    eig(full(Prr1))
    
    Xn1=[state1(1:6);params1;state1(7:14);params2];
    
    
    %%
    % estimated state:
%     y0=yn;
    X0=Xn1;
    x0=xn;
    
    feats_HK0=feats_HKn;
    
    Prr0=0.5*(Prr1+Prr1');
    
    frameLeftGray0=frameLeftGrayn;
    filt_match10=filt_match1n;
    filt_match20=filt_match2n;
    
    live_processing;
    
    %% message information:
    fprintf('\t Iteration nr. %d  \n', loop);
    %     fprintf('------------------------------------------\n');
    %     fprintf('->\t Simulation step          : %.3f  [s]\n', t_sim);
    %     fprintf('->\t Filter propagation step  : %.3f  [s]\n', t_prop);
    %     fprintf('->\t Matching measure step    : %.3f  [s]\n', t_match);
    %     fprintf('->\t Filter update step       : %.3f  [s]\n', t_update);
    %     fprintf('->\t Landmarks Initialization : %.3f  [s]\n', t_init);
    %     fprintf('->\t Iteration storage        : %.3f  [s]\n', t_storage);
    %     fprintf('------------------------------------------\n');
    %     fprintf('\t Tot.                     : %.3f  [s]\n', t_storage+t_sim+t_prop+t_match+t_update+t_init);
    
    
end
