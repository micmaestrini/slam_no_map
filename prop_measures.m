function [y,lmkinfo]=prop_measures(Xn,S0,Prr,Prm,Pmm,cam_params,qc,Rn,lmkinfo)
    % Function that estimates measures available at the state x^-_{k+1}.
    % Inputs: 
    % X0    : state of the filter after previous correction (i.e. x^+_k);
    % S0    : current knowledge of landmarks [nm x 3];
    % Hf    : sensor horizontal size[m];
    % Vf    : sensor vertical size[m];
    % foc   : focal length [m];
    % b     : baseline between sensors [m];
    % 
    % Outputs:
    % Yn    : estimated measures (including out of FOV) [4 x nm];
    % vis   : vector of indexes of visible landmarks (subset of current map);

    % extraction of position and relative attitude from full state:
        x0=Xn(1);
        y0=Xn(2);
        z0=Xn(3);
        f0=Xn(9);
        s=Xn(14:16);
        s1=s(1);
        s2=s(2);
        s3=s(3);

    % DCM that rotates from LVLH to chaser body frame definition:
        C_BI=quat2dcm(qc);
        C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL=C_BI*C_LI';
        qc=dcm2quat(C_BL);
        sc=qc(2:end)/(1+qc(1));

    % DCM that rotates from target to chaser frame (D'):
        skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
        D=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;
        D=D';

    % rotation of current landmarks in the chaser reference frame:
        P_i=D*S0'+C_BL*[x0;y0;z0];
        xi=P_i(1,:);
        yi=P_i(2,:);
        zi=P_i(3,:);

    % stereo measurement equation:
        h=[cam_params.u0-cam_params.alpha_u*xi./zi;cam_params.v0-cam_params.alpha_v*yi./zi;cam_params.alpha_u*cam_params.b./zi];
%         scatter(h(1,:),h(2,:))

    % visibility is assumed if projection is expected inside FOV for both
    % cameras:
        vis=find(abs(h(1,:))<cam_params.hpix/2 & abs(h(2,:))<cam_params.vpix/2);
        vis=linspace(1,size(h,2),size(h,2));
        lmkinfo.counter_prop=lmkinfo.counter_prop+1;
        
        y.h=reshape(h,[],1);
        y.vis=vis;
        n=length(vis);
        y.n=n;
        
        
        % n is the size of visible estimated landmarks:
        
    % nh is the size of one landmark measurement:
        nh=3;
    % extraction of attitude and position from current state:

    % initialization of dummy vectors that store the values of the rows,
    % columns and content of the sparse partial of the measurement.
    % 36 is obtained by counting the number of outputs of function generated
    % with symbolic:
        ROWS=zeros(n*(18+9),1);
        COLUMNS=zeros(n*(18+9),1);
        VALUES=zeros(n*(18+9),1);
        y.sigma=zeros(n,1);

    % loop over all visible points:
        for i=1:n
                h_pi=H_Pi(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,S0(vis(i),1),S0(vis(i),2),S0(vis(i),3),s1,s2,s3,sc(1),sc(2),sc(3),x0,y0,z0);
                h_xi=H_x(cam_params.alpha_u,cam_params.alpha_v,cam_params.b,S0(vis(i),1),S0(vis(i),2),S0(vis(i),3),s1,s2,s3,sc(1),sc(2),sc(3),x0,y0,z0);
                h_i=[h_xi,h_pi];
                
                columns_extraction_i=reshape(3*vis(i)-2+(0:2)',[],1);
                Pn_i=[Prr,Prm(:,columns_extraction_i);Prm(:,columns_extraction_i)',Pmm(columns_extraction_i,columns_extraction_i)];

                y.sigma(i)=det(h_i*Pn_i*h_i'+Rn);
                
            % unpacking sparse jacobian in rows, columns and values:
                [r1,c1,v1]=find(h_xi);
                [r2,c2,v2]=find(h_pi);
            % rows are incremented to stack all matrix in one:
                r1=nh*(i-1)+r1;
                r2=nh*(i-1)+r2;
            % columns are corrected to account for fixed number of state:
                c2=14+3*(i-1)+c2;
            % storage of values of one iteration;
%             length([r1;r2])
                ROWS(27*i-26:27*i)=[r1;r2];
                COLUMNS(27*i-26:27*i)=[c1;c2];
                VALUES(27*i-26:27*i)=[v1;v2];
        end

    % assembly of full scale H_ matrix with assumed output size:
        y.H=sparse(ROWS,COLUMNS,VALUES,nh*n,14+3*n);
        
    % assembly of full scale Pn matrix excluding points outside FOV:
        columns_extraction=reshape(3*vis-2+(0:2)',[],1);
        Pn=[Prr,Prm(:,columns_extraction);Prm(:,columns_extraction)',Pmm(columns_extraction,columns_extraction)];

    % assembly of full scale measurement noise matrix:
        RnCell = repmat({Rn}, 1, n);
    % it is a blockdiagonal as each measure is independent from another:
        BigR = blkdiag(RnCell{:});
        if isempty(BigR)
            BigR=zeros(0,0);
        end
    % Computation of innovation covariance: 
        y.HPH=y.H*Pn*y.H'+BigR;
end