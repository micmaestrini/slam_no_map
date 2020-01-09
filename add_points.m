function [S_new,Prm_new,Pmm_new]=add_points(X0,y0,Prr0,Prm0,Pmm0,R,cam_params,indexPairs)
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
        new=double(indexPairs(:,1)');
        indexes=reshape((new-1)*3+[1:3]',[],1);
        unmatched_measures=y0.z(indexes);
    % n is the size of unmatched measures:
        n=length(new);
        %         unmatched_measures(:,match(:,2))=[];

    % inverse measure function:
        f0=X0(9);
        x0=X0(1);
        y0=X0(2);
        z0=X0(3);
        s0=X0(14:16);

        sc=X0(22:24);
        skew_sc=[0,-sc(3),sc(2);sc(3),0,-sc(1);-sc(2),sc(1),0];
        C_BI=eye(3)+8*(skew_sc*skew_sc)/(1+transpose(sc)*sc)^2-4*(1-transpose(sc)*sc)/(1+transpose(sc)*sc)^2*skew_sc;
        C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL=C_BI*C_LI';
        qc=dcm2quat(C_BL);
        sc=qc(2:end)/(1+qc(1));





        
        P_i =@(m) G_fun(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,m(3),s0(1),s0(2),s0(3),sc(1),sc(2),sc(3),m(1),cam_params.u0,m(2),cam_params.v0,x0,y0,z0);
        dg_pi = @(m) G_Pi(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,m(3),s0(1),s0(2),s0(3),m(1),cam_params.u0,m(2),cam_params.v0);
        dg = @(m) G_x(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,m(3),s0(1),s0(2),s0(3),sc(1),sc(2),sc(3),m(1),cam_params.u0,m(2),cam_params.v0,x0,y0,z0);

    % initialization of measured points:
        Pi=zeros(3,n);

    % ns is the size of state output of inverse measure function:
        ns=3;
        nh=3;
    % initialization of cell structure to contain all jacobians wrt measure
    % (will be block diagonal):
        GPi=cell(n,1);
    % initialization of measured points:

    % dummy variable so store rows, columns and values of the entrie jacobian
    % that needs to be assembled as a sparse matrix:
        ROWS=zeros(n*(18),1);
        COLUMNS=zeros(n*(18),1);
        VALUES=zeros(n*(18),1);
    % loop over unmatched measures:
        for i=1:n
                rows=(i-1)*nh+[1:nh]';

            % add new landmark in 3d:
                pi_new = P_i(unmatched_measures(rows));
            % compute jacobian matrix wrt state:
                gri = dg(unmatched_measures(rows));
            % compute jacobian matrix wrt measure:
                gpi = dg_pi(unmatched_measures(rows));

            % store measure jacobian in cell array:
                GPi{i}=gpi;
            % store new point in array:
                Pi(:,i)=pi_new;

            % unpack state jacobian in rows, columns and values.
                [r1,c1,v1]=find(gri);

            % stack values remembering to increase the rows ad each iteration:
                ROWS(18*i-17:18*i)=ns*(i-1)+r1;
                COLUMNS(18*i-17:18*i)=c1;
                VALUES(18*i-17:18*i)=v1;
        end

    % assemble complete measure hessian with suggested otuput size:
        GR=sparse(ROWS,COLUMNS,VALUES,3*n,14);

    % assemble jacobian of measure as block diagonal:
        GP=blkdiag(GPi{:});

    % assemble of block diagonal matrix for measurement noise:
        RnCell = repmat({R}, 1, n);
        BigR = blkdiag(RnCell{:});

    % creation of new blocks for covariance matrix. Pll is new landmarks only.
    % Plr is new landmarks with state and Plm is new landmark with all other
    % landmarks:
        Pll=GR*Prr0*GR'+GP*BigR*GP';
        Plr=GR*Prr0;
        Plm=GR*Prm0;

    % assembly of useful terms only:
        Prm_new=[Prm0,Plr'];
        Pmm_new=[Pmm0,Plm';Plm,Pll];
        S_new=Pi';

end