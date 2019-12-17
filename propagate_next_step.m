function [Xn,Yn,Prrn,Prmn,Pmmn,lmkinfo]=propagate_next_step(X0,Prr0,Prm0,Pmm0,Q0,dt,S0,cam_params,qc,R,lmkinfo)
    % Function that propagates from the state x^+_k to the state x^-_{k+1}.
    % Inputs: 
    % X0    : state of the filter after previous correction (i.e. x^+_k);
    % Prr0  : covariance submatrix for state only terms at k_+;
    % Prm0  : covariance submatrix for mixed state-landmarks terms at k_+;
    % Pmm0  : covariance submatrix for  landmarks only terms at k_+;
    % Q0    : Process noise matrix during k_ transition;
    % dt    : time interval of propagation [s];
    % S0    : current knowledge of landmarks [nm x 3];
    % Hf    : sensor horizontal size[m];
    % Vf    : sensor vertical size[m];
    % foc   : focal length [m];
    % b     : baseline between sensors [m];
    % 
    % Outputs:
    % Xn    : state of the filter before correction (i.e. x^-_{k+1});
    % Yn    : estimated measures (including out of FOV) [4 x nm];
    % Prrn  : covariance submatrix for state only terms at k+1_-;
    % Prmn  : covariance submatrix for mixed state-landmarks terms at k+1_-;
    % Pmmn  : covariance submatrix for  landmarks only terms at k+1_-;
    % vis   : vector of indexes of visible landmarks (subset of current map);

    %% propagate next state:
    % step to propagate the states inside the filter:
        [Xn]=prop_states(X0,dt);

    %% propagate covariance:
    % covariance propagation step:
        [Prrn,Prmn,Pmmn]=prop_covariance(X0,dt,Prr0,Prm0,Pmm0,Q0);

    %% estimate new measures at estimated new state:
    % estimate of measures and visibility (inside FOV):
        [Yn,lmkinfo]=prop_measures(Xn,S0,Prrn,Prmn,Pmmn,cam_params,qc,R,lmkinfo);

end