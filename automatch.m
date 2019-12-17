function [match,Hn,Zn]=automatch(index_meas,vis,map_visibility,Xn,S0,Prr,Prm,Pmm,Rn,b,foc,qc)
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
    

    index_of_sim_points=map_visibility(vis);
    right_measures=find(ismember(index_of_sim_points,index_meas));
    meas_index_matched=find(ismember(index_meas,index_of_sim_points(right_measures)));
    if isempty(meas_index_matched)
        meas_index_matched=zeros(0,1);
    end
    if isempty(right_measures)
        right_measures=zeros(0,1);
    end
    
    match=[vis(right_measures),meas_index_matched];
   
    % n is the size of visible estimated landmarks:
        n=length(vis);
    % nh is the size of one landmark measurement:
        nh=4;
    % extraction of attitude and position from current state:
        s0=Xn(14:16);
        x0=Xn(1);
        y0=Xn(2);
        z0=Xn(3);
        f0=Xn(9);
        
        
        C_BI=quat2dcm(qc);
        C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL=C_BI*C_LI';
        qc=dcm2quat(C_BL);

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

end



