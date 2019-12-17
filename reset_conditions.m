function [x0,X0,Prr0,Prm0,Pmm0,S0]=reset_conditions(xn,Xn,Prrn,Prmn,Pmmn,S_new,Sn)
    % Function that resets filter states and covariances to continue iterating.
    % Inputs:
    % xn     : simulator state at k+1 [4 x nmeas];
    % Xn     : state of the filter after correction (i.e. x^+_k+1);
    % Prrn   : covariance submatrix for state only terms at k+1_+;
    % Prmn   : covariance submatrix for mixed state-landmarks terms at k+1_+;
    % Pmmn   : covariance submatrix for  landmarks only terms at k+1_+;
    % S_new  : newly initialized landmarks [n_new,3];
    % Sn     : landmarks estimate at k+1_+;
    %
    % Outputs:
    % x0     : new simulator state [4 x nmeas];
    % X0     : new filter state;
    % Prr0   : new covariance submatrix for state only;
    % Prm0   : new covariance submatrix for mixed state-landmarks;
    % Pmm0   : new covariance submatrix for  landmarks only;
    % S0     : new landmarks map;

    % simulator state:
        x0=xn;
    % estimated state:
        X0=Xn;
    % estimated new covariance:
        P=[Prrn,Prmn;Prmn',Pmmn];
        P=0.5*(P+P');
        Prr0=P(1:14,1:14);
        Prm0=P(1:14,15:end);
        Pmm0=P(15:end,15:end);
        
    % estimate new set of points:
        S0=[Sn;S_new'];
end