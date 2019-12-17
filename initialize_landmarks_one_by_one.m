function [Pi,Prm,Pmm]=initialize_landmarks_one_by_one(X,yn,Prr,Prm,Pmm,R,cam_params,qc)
    % Function that adds unmatched landmarks measures to the state map.
    % Inputs: 
    % X      : state of the filter after correction (i.e. x^+_k+1);
    % yn     : simulated measures [4 x nmeas];
    % match  : vector of size [nmatch,2], whose columns contain indexes of
    % estimated and real measures that were matched;
    % b      :  baseline between cameras along y [m];
    % foc    :  focal length for camera [m];
    % Prr    : covariance submatrix for state only terms at k+1_+;
    % Prm    : covariance submatrix for mixed state-landmarks terms at k+1_+;
    % Pmm    : covariance submatrix for  landmarks only terms at k+1_+;
    % Rn     : measurement noise matrix of size 4x4;
    % 
    % Outputs:
    % Pi     : newly initialized landmarks [n_new,3];
    % Prm    : covariance submatrix for mixed state-landmarks after addition;
    % Pmm    : covariance submatrix for  landmarks only after addition;


    % initialization of unmatched measures by subtraction from matched
    % measurement:
        unmatched_measures=yn;

    % inverse measure function:
        f0=X(9);
        x0=X(1);
        y0=X(2);
        z0=X(3);
        s0=X(14:16);
        
        C_BI=quat2dcm(qc);
        C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL=C_BI*C_LI';
        qc=dcm2quat(C_BL);
        sc=qc(2:end)/(1+qc(1));
        
        P_i =@(m) G_fun(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,m(3),s0(1),s0(2),s0(3),sc(1),sc(2),sc(3),m(1),cam_params.u0,m(2),cam_params.v0,x0,y0,z0);
        dg_pi = @(m) G_Pi(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,m(3),s0(1),s0(2),s0(3),m(1),cam_params.u0,m(2),cam_params.v0);
        dg =G_x(s0(1),s0(2),s0(3),sc(1),sc(2),sc(3));
        
    % n is the size of unmatched measures:
        n=min(5,size(unmatched_measures,2));

    % initialization of measured points:
        Pi=zeros(3,n);

    % loop over unmatched measures:
        for i=1:n
            % compute jacobian matrix wrt measure:
                gpi = dg_pi(unmatched_measures(:,i));

            % store new point in array:
                Pi(:,i)=P_i(unmatched_measures(:,i));
                
                Plli=dg*Prr*dg'+gpi*R*gpi';
                Plri=dg*Prr;
                Plmi=dg*Prm;
                Prm=[Prm,Plri'];
                Pmm=[Pmm,Plmi';Plmi,Plli];
        end

end