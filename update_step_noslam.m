function [Prrn,Prmn,Pmmn,Xn,Sn]=update_step_noslam(Xn,Prr,Prm,Pmm,yn,Yn,S0)
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
    zk = yn.z(1:15);
    hk = Yn.h;
    % Rk = observations.R(indi,indi);
    Ck = Yn.HPH;
    Hk = Yn.H;
    
    % remove r, dr, theta,dtheta from state and save them in dummy variable:
        state=Xn(1:6);
        params=Xn(7:18);
        
    % assemble full covariance matrix:
        P0=[Prr,Prm;Prm',Pmm];

    % define augmeted state:
        X=[state;reshape(S0',[],1)];

    %% Kalman update:
    % Kallman gain computation:
        K=P0*Hk'/Ck;
        
    % state update:
        state_n=X+K*(zk-hk);
        d2 = mahalanobis (zk-hk, Ck)
    % covariance update:
        Pn=P0-K*Ck*K';
        
    %% output reshaping:
    % updated state reaugmented with parameters:
        Xn=[state_n(1:6);params];
        state_n(1:6)-state
    % reshaping landmarks after update:
        Sn=reshape(state_n(7:end),3,[])';

    % extraction of relevant submatrices:
        Prrn=Pn(1:6,1:6);
        Prmn=Pn(1:6,7:end);
        Pmmn=Pn(7:end,7:end);

end