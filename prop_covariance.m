function [Prr1,Prm1,Pmm1]=prop_covariance(X0,dt,Prr0,Prm0,Pmm0,Q0)
    % This function propagates the states and landmark covariance :
    % Inputs:
    % X0    : initial state of filter before propagation k_+;
    % Prr0  : covariance submatrix for state only terms at k_+;
    % Prm0  : covariance submatrix for mixed state-landmarks terms at k_+;
    % Pmm0  : covariance submatrix for  landmarks only terms at k_+;
    % Q0    : Process noise matrix during k_ transition;
    %
    % Outputs:
    % Prr1  : covariance submatrix for state only terms at k+1_-;
    % Prm1  : covariance submatrix for mixed state-landmarks terms at k+1_-;
    % Pmm1  : covariance submatrix for  landmarks only terms at k+1_-;

    % state transition matrix computation:
        Phi0=STM(X0,dt);
    % state only covariance propagation:
        Prr1=Phi0*Prr0*Phi0'+Q0;
    % mixed state landmarks covariance propagation:
        Prm1=Phi0*Prm0;
    % the covariance propagation is independent of the collected points and
    % does not influence their pure covariance Pmm, just the mixed state terms.
        Pmm1=Pmm0;

end