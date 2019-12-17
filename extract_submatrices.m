function [P_,H_,Z_,dz_]=extract_submatrices(match,vis,Zn,Prr,Prm,Pmm,Hn,yn,Yn)
    % Function that extracts submatrices of the innovation and measure jacobian
    % from their full form in order to have the correct size for update.
    % Moreover extracts also the state covariance required and the innovation
    % vector.
    % Inputs:
    % match  : vector of size [nmatch,2], whose columns contain indexes of
    % estimated and real measures that were matched;
    % vis    : vector of indexes of visible landmarks (subset of current map);
    % Zn     : innovation covariance matrix;
    % Prr    : covariance submatrix for state only terms at k+1_-;
    % Prm    : covariance submatrix for mixed state-landmarks terms at k+1_-;
    % Pmm    : covariance submatrix for  landmarks only terms at k+1_-;
    % Hn     : measurement jacobian;
    % yn     : simulated measures [4 x nmeas];
    % Yn     : estimated measures (including out of FOV) [4 x nm];
    % S0     : current knowledge of landmarks [nm x 3];
    %
    % Outputs:
    % P_     : covariance submatrix of correct size for update;
    % H_     : jacobian submatrix of measurements of correct size for update;
    % Z_     : covariance submatrix for innovation of correct size for update;
    % dz_    : innovation of correct size for update;

    % Hn and Zn are already based on visible points only:
    % Zn size is (4xvis)x(4xvis)
    % finds indexes of matched points in the visible vector:
        visibility=find(ismember(vis,match(:,1)));
        if isempty(visibility)
            visibility=zeros(0,1);
        end
    % size of measure vector for one point:
        nh=4;
    % output of the inverse measurement equation for one point:
        ns=3;

    % extract index of rows to be kept from H and Z matrices which do not have
    % m size, but vis size (computed from the visible estimate points):
        keep_vis_index_meas=reshape([nh*visibility-(nh-1)+[0:(nh-1)]]',[],1);
        keep_vis_index_state=reshape([ns*visibility-(ns-1)+[0:(ns-1)]]',[],1);
        keep_vis_index_covmat=reshape([ns*match(:,1)-(ns-1)+[0:(ns-1)]]',[],1);
    % extract submatrices using generated indexes:
        Z_=Zn(keep_vis_index_meas,keep_vis_index_meas);
    % notice shift of 14 because H has always the state components in the first
    % 14 columns:
        H_=Hn(keep_vis_index_meas,[(1:6)';6+keep_vis_index_state]);
    % innovation computation based on match:
        dz_=reshape(yn(:,match(:,2))-Yn(:,match(:,1)),[],1);
    % covariance matrix truncated to have correct size:
        P_=[Prr,Prm(:,keep_vis_index_covmat);Prm',Pmm(:,keep_vis_index_covmat)];

end