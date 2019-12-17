%% partials chaser angular vel

clear
clc
close all

syms f vf b0 b1 b2 b3 db0 db1 db2 db3 ddb0 ddb1 ddb2 ddb3 ddf

C11=1-2*b2^2-2*b3^2;
C12=2*(b1*b2+b0*b3);
C13=2*(b1*b3-b0*b2);
C21=2*(b1*b2-b0*b3);
C31=2*(b1*b3+b0*b2);
C22=1-2*b1^2-2*b3^2;
C23=2*(b2*b3+b0*b1);
C32=2*(b2*b3-b0*b1);
C33=1-2*b1^2-2*b2^2;
C_BL=[C11,C12,C13;C21,C22,C23;C31,C32,C33];
% C_BL=[C11,C21,C31;C12,C22,C32;C13,C23,C33];
% C_IL:
C_LI=[cos(f),sin(f),0;-sin(f),cos(f),0;0,0,1];
% C_BI
C_BI=simplify(C_BL*C_LI);
C_IB=transpose(C_BI);


C_IB0=diff(C_IB,b0);
C_IB1=diff(C_IB,b1);
C_IB2=diff(C_IB,b2);
C_IB3=diff(C_IB,b3);
C_IBteta=diff(C_IB,f);


C_t=simplify(C_IB0*db0+C_IB1*db1+C_IB2*db2+C_IB3*db3+C_IBteta*vf);

W=-C_t*C_IB;
omega=([-W(2,3);W(1,3);-W(1,2)]);


omega_0=simplify(diff(omega,b0));
omega_1=simplify(diff(omega,b1));
omega_2=simplify(diff(omega,b2));
omega_3=simplify(diff(omega,b3));
omega_d0=simplify(diff(omega,db0));
omega_d1=simplify(diff(omega,db1));
omega_d2=simplify(diff(omega,db2));
omega_d3=simplify(diff(omega,db3));
omega_f=simplify(diff(omega,f));
omega_vf=simplify(diff(omega,vf));

omega_t=omega_0*db0+omega_1*db1+omega_2*db2+omega_3*db3+...
        omega_d0*ddb0+omega_d1*ddb1+omega_d2*ddb2+omega_d3*ddb3+...
        omega_f*vf+omega_vf*ddf;


omega_b=jacobian(omega,[b0;b1;b2;b3]);
omega_db=jacobian(omega,[db0;db1;db2;db3]);
domega_b=jacobian(omega_t,[b0,b1,b2,b3]);
domega_db=jacobian(omega_t,[db0,db1,db2,db3]);
domega_ddb=jacobian(omega_t,[ddb0,ddb1,ddb2,ddb3]);


% matlabFunction(omega,'File','wc_fun','Sparse',true);
% matlabFunction(omega_t,'File','dwc_fun','Sparse',true);
% 
% matlabFunction(omega_b,'File','STM_wc_b','Sparse',true);
% matlabFunction(omega_db,'File','STM_wc_db','Sparse',true);
% 
% matlabFunction(domega_b,'File','STM_dwc_b','Sparse',true);
% matlabFunction(domega_db,'File','STM_dwc_db','Sparse',true);
% matlabFunction(domega_ddb,'File','STM_dwc_ddb','Sparse',true);





