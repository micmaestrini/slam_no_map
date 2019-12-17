function [Prrn,Prmn,Pmmn,Xn,Sn]=correction_one_by_one(b,foc,Xn,Prr,Prm,Pmm,yn,Yn,S0,R,vis,match,qc)
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
    
    

        f0=Xn(9);
        s0=Xn(14:16);
        
        C_BI=quat2dcm(qc);
        C_LI=[cos(f0),sin(f0),0;-sin(f0),cos(f0),0;0,0,1];
        C_BL=C_BI*C_LI';
        qc=dcm2quat(C_BL);
    
    % output of the inverse measurement equation for one point:
        ns=3;
        
        P=[Prr,Prm;Prm',Pmm];
        X=[Xn(1:6);reshape(S0',[],1)];
    
    for i=1:size(match,1)
        
        x0=X(1);
        y0=X(2);
        z0=X(3);
        Prr=P(1:6,1:6);
        Prm=P(1:6,7:end);
        Pmm=P(7:end,7:end);
        
        dzi=yn(:,match(i,2))-Yn(:,match(i,1));
    % partial of the measure wrt the state:
        h_xi=H(b,foc,S0(vis(i),1),S0(vis(i),2),S0(vis(i),3),qc(1),qc(2),qc(3),qc(4),s0(1),s0(2),s0(3),x0,y0,z0);
    % partial of the measure wrt the landmark position:
        h_pi=H_pi(b,foc,S0(vis(i),1),S0(vis(i),2),S0(vis(i),3),qc(1),qc(2),qc(3),qc(4),s0(1),s0(2),s0(3),x0,y0,z0);
        keep_vis_index_covmat=[ns*match(i,1)-(ns-1)+[0:(ns-1)]]';
          
        Hi=[h_xi,h_pi];
        Pi=[Prr,Prm(:,keep_vis_index_covmat);Prm(:,keep_vis_index_covmat)',Pmm(keep_vis_index_covmat,keep_vis_index_covmat)];
        Zi=Hi*Pi*Hi'+R;
%         Zi-Zi'
        if (dzi'*(Zi\dzi)<0)
           disp('ciao############################################');
        end
        Ki=[Prr,Prm(:,keep_vis_index_covmat);Prm',Pmm(:,keep_vis_index_covmat)]*Hi'*inv(Zi);
        
        if (dzi'*(Zi\dzi))<9       
            X=X+Ki*dzi;
            P=P-Ki*Zi*Ki';
            S0=reshape(X(7:end),3,[])';
        end
               
    end
    
    Xn(1:6)=X(1:6);
    Sn=reshape(X(7:end),3,[])';
    Prrn=P(1:6,1:6);
    Prmn=P(1:6,7:end);
    Pmmn=P(7:end,7:end);

end