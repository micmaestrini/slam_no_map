function [simulator_history,state_history,landmarks_map,covariance_mat,visibility_check]=Store_Iteration(xn,yn,...
    Xk,Sk,Yk,Prrk,Prmk,Pmmk,...
    Prrk1,Prmk1,Pmmk1,Xk1,Sk1,...
    Prm_new,Pmm_new,S_new,...
    simulator_history,state_history,landmarks_map,covariance_mat,counter,visibility_check)

% tic;

%% Simulator Processing:
simulator_history.state(:,counter)=xn;
simulator_history.measures{counter}=yn;

%% State Processing:
state_history.xk(:,counter)=Xk;
state_history.xk1(:,counter)=Xk1;

%% Filter Measures Processing:
visibility_check.filter_meas{counter}=Yk;

%% Landmarks Processing:
landmarks_map.Sk{counter}=Sk;
landmarks_map.Sk1{counter}=Sk1;
landmarks_map.new{counter}=S_new;

%% Covariance Processing:
covariance_mat.Prrk{counter}=Prrk;
covariance_mat.Prmk{counter}=Prmk;
covariance_mat.Pmmk{counter}=Pmmk;

covariance_mat.Prrk1{counter}=Prrk1;
covariance_mat.Prmk1{counter}=Prmk1;
covariance_mat.Pmmk1{counter}=Pmmk1;

covariance_mat.Prm_augmented{counter}=Prm_new;
covariance_mat.Pmm_augmented{counter}=Pmm_new;

% t_f=toc;

%% Times Processing:
% T=[t_sim,t_match,t_update,t_init,t_f];
% times_history.iter(counter,:)=T;


end