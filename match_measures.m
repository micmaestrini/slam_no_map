function [match,Hn,Zn]=match_measures(Xn,yn,Yn,Prr,Prm,Pmm,Rn,S0,vis,b,foc,qc)
    % Function that matches measurements obtained from simulator with expected
    % measurements from filter propagation.
    % Inputs: 
    % Xn    : state of the filter before correction (i.e. x^-_k+1);
    % yn    : simulated measures [4 x nmeas];
    % Yn    : estimated measures (including out of FOV) [4 x nm];
    % Prr   : covariance submatrix for state only terms at k+1_-;
    % Prm   : covariance submatrix for mixed state-landmarks terms at k+1_-;
    % Pmm   : covariance submatrix for  landmarks only terms at k+1_-;
    % Rn    : measurement noise matrix of size 4x4;
    % S0    : current knowledge of landmarks [nm x 3];
    % vis   : vector of indexes of visible landmarks (subset of current map);
    % b     : baseline between sensors [m];
    % foc   : focal length [m];
    % 
    % Outputs:
    % Hn    : measurement jacobian;
    % Zn    : innovation covariance matrix;
    % match : vector of size [nmatch,2], whose columns contain indexes of
    % estimated and real measures that were matched.

    % n is the size of visible estimated landmarks:
        n=length(vis);
    % nh is the size of one landmark measurement:
        nh=4;
    % extraction of attitude and position from current state:
        s0=Xn(14:16);
        x0=Xn(1);
        y0=Xn(2);
        z0=Xn(3);

    % initialization of dummy vectors that store the values of the rows,
    % columns and content of the sparse partial of the measurement.
    % 36 is obtained by counting the number of outputs of function generated
    % with symbolic:
        ROWS=zeros(n*(12+12),1);
        COLUMNS=zeros(n*(12+12),1);
        VALUES=zeros(n*(12+12),1);
    % loop over all visible points:
        for i=1:n
            % partial of the measure wrt the state:
                h_xi=H(b,foc,S0(vis(i),1),S0(vis(i),2),S0(vis(i),3),qc(1),qc(2),qc(3),qc(4),s0(1),s0(2),s0(3),x0,y0,z0);
            % partial of the measure wrt the landmark position:
                h_pi=H_pi(b,foc,S0(vis(i),1),S0(vis(i),2),S0(vis(i),3),qc(1),qc(2),qc(3),qc(4),s0(1),s0(2),s0(3),x0,y0,z0);
            % unpacking sparse jacobian in rows, columns and values:
                [r1,c1,v1]=find(h_xi);
                [r2,c2,v2]=find(h_pi);
            % rows are incremented to stack all matrix in one:
                r1=nh*(i-1)+r1;
                r2=nh*(i-1)+r2;
            % columns are corrected to account for fixed number of state:
                c2=6+3*(i-1)+c2;
            % storage of values of one iteration;
                ROWS(24*i-23:24*i)=[r1;r2];
                COLUMNS(24*i-23:24*i)=[c1;c2];
                VALUES(24*i-23:24*i)=[v1;v2];
        end

    % assembly of full scale H_ matrix with assumed output size:
        Hn=sparse(ROWS,COLUMNS,VALUES,nh*n,6+3*n);

    % assembly of full scale Pn matrix excluding points outside FOV:
        columns_extraction=reshape(3*vis'-2+(0:2)',[],1);
        Pn=[Prr,Prm(:,columns_extraction);Prm(:,columns_extraction)',Pmm(columns_extraction,columns_extraction)];

    % assembly of full scale measurement noise matrix:
        RnCell = repmat({Rn}, 1, n);
    % it is a blockdiagonal as each measure is independent from another:
        BigR = blkdiag(RnCell{:});

    % Computation of innovation covariance: 
        Zn=Hn*Pn*Hn'+BigR;

    % initialization of matching array:
        match=zeros(0,2);
    % eval set is the array containing index of visible points wrt array S0.
        eval_set=1:size(yn,2);

    % loop over estimated visible points to try to match them to the measures:
        for i=1:n
            % extract submatrix of the innovation that characterize one measure (2
            % images):
                Pr=Zn(4*i-3:4*i,4*i-3:4*i);
                if issparse(Pr)
                    Pr=full(Pr);
                end
                    
                [R,D]=svd(Pr);
                invPr=R*diag(1./diag(D))*R';

            % Squared mahalanobis distance function:
                F=@(dz) dot(dz,invPr*dz,1);
            % innovation from one visible estimated measure to all the available
            % measurements:
                dz=yn(:,eval_set)-Yn(:,vis(i));
            % squared mahalanobis distance of each innovation:
                variance=F(dz);
            % find indexes with smaller than 3sigma variance:
                matched_ind=find(variance<1);

            % if there are more than one points inside the uncertainty ellipsoid of
            % 3 sigma about one measure:
                if length(matched_ind)>1
                    % extract the position of the minimum value of the variance of all
                    % matched points:
                        [~,pos]=min(variance(matched_ind));
                    % the index of this minimum is usewd to extract the value inside
                    % matched_ind, which corresponds to the index on the eval set. The
                    % index on the eval set is hence extracted.
                        ind_match=eval_set(matched_ind(pos));
                    % the size of the match vector is incremented by the new couple
                    % (estimate, measure):
                        match=[match;vis(i), ind_match];
                    % the measure is removed from the evaluation set:
                        eval_set(matched_ind(pos))=[];

                % if there is only one point inside the uncertainty ellipsoid of 3 
                % sigma about one measure:
                elseif length(matched_ind)==1
                    % the index of the matching is used to extract the index of the
                    % measure from the evaluation set:
                        ind_match=eval_set(matched_ind);
                    % the size of the match vector is incremented by the new couple
                    % (estimate, measure):
                        match=[match;vis(i), ind_match];
                    % the measure is removed from the evaluation set:
                        eval_set(matched_ind)=[];
                end

        end

end



