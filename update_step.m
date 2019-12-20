function [Prrn,Prmn,Pmmn,Xn,Sn,lmkinfo,feats_list]=update_step(Xn,Prr,Prm,Pmm,yn,Yn,S0,H,Opt,lmkinfo,feats_list)
    % Function that updates states and covariances from step k+1_- to k+1_+.
    % Inputs: 
    % X0     : state of the filter before correction (i.e. x^-_k+1);
    % Prr    : covariance submatrix for state only terms at k+1_-;
    % Prm    : covariance submatrix for mixed state-landmarks terms at k+1_-;
    % Pmm    : covariance submatrix for  landmarks only terms at k+1_-;
    % Zn     : innovation covariance matrix;
    % Hn     : measurement jacobian;
    % yn     : simulated measures [4 x nmeas];
    % Yn     : estimated measures (including out of FOV) [4 x nm];
    % S0     : current knowledge of landmarks [nm x 3];
    % vis    : vector of indexes of visible landmarks (subset of current map);
    % match  : vector of size [nmatch,2], whose columns contain indexes of
    % estimated and real measures that were matched.
    % 
    % Outputs:
    % Prrn   : covariance submatrix for state only terms at k+1_+;
    % Prmn   : covariance submatrix for mixed state-landmarks terms at k+1_+;
    % Pmmn   : covariance submatrix for  landmarks only terms at k+1_+;
    % Xn     : state of the filter after correction (i.e. x^+_k+1);
    % Sn     : updated knowledge of landmarks [nm x 3];

    % extraction of useful submatrices of correct size:
    % i meas, j prop
    [~,i, j] = find(H);
    
    indexMatched=matchFeatures(feats_list(j,:),yn.feats(i,:));
    keep_meas=indexMatched(:,1)==indexMatched(:,2);
    i=i(keep_meas);
    j=j(keep_meas);
    
    feats_list(j,:)=yn.feats(i,:);    
    
    
    lmkinfo.counter_meas(j)=lmkinfo.counter_meas(j)+1;
    
    [~,pos]= sort(Yn.sigma(j),'desc');
%     pos=randperm(length(j));
%     i(pos)
%     j(pos)

    indi=reshape(3*i(pos)-2+[0:2]',[],1);
    indj=reshape(3*j(pos)-2+[0:2]',[],1);
    
    indi=indi(1:min(Opt.Nup*3,length(indi)));
    indj=indj(1:min(Opt.Nup*3,length(indj)));

    zk = yn.z(indi);
    hk = Yn.h(indj);
    % Rk = observations.R(indi,indi);
    Ck = Yn.HPH(indj,indj);
    Hk = Yn.H(indj,:);
    
    % remove r, dr, theta,dtheta from state and save them in dummy variable:
        state=Xn([1:6,11:18]);
        params1=Xn(7:10);
        params2=Xn(19:24);
        
    % assemble full covariance matrix:
        P0=[Prr,Prm;Prm',Pmm];

    % define augmeted state:
        X=[state;reshape(S0',[],1)];

    %% Kalman update:
    % Kallman gain computation:
        K=P0*Hk'/Ck;
    % state update:
        state_n=X+K*(zk-hk);
    % covariance update:
        Pn=P0-K*Ck*K';
        
    %% output reshaping:
    % updated state reaugmented with parameters:
        Xn=[state_n(1:6);params1;state_n(7:14);params2];
    % reshaping landmarks after update:
        Sn=reshape(state_n(15:end),3,[])';

    % extraction of relevant submatrices:
        Prrn=Pn(1:14,1:14);
        Prmn=Pn(1:14,15:end);
        Pmmn=Pn(15:end,15:end);

end