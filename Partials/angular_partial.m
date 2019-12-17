%% partials angular vel

clear
clc
close all

syms s1 s2 s3 k1 k2 wx wy wz wcx wcy wcz dwcx dwcy dwcz


J=diag([exp(k1),1,exp(k2)]);
invJ=diag([1/exp(k1),1,1/exp(k2)]);
% D da C a T!

s=[s1;s2;s3];
skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
D=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;

w=[wx;wy;wz];
wc=[wcx;wcy;wcz];
dwc=[dwcx;dwcy;dwcz];


dw=invJ*(-cross(w,J*w)+J*cross(w,D*wc)-J*D*dwc-cross(D*wc,J*D*wc)-cross(w,J*D*wc)-cross(D*wc,J*w));

dw=simplify(dw);

partials_w=simplify(jacobian(dw,[wx;wy;wz]));
partials_s=simplify(jacobian(dw,[s1;s2;s3]));
partials_wc=simplify(jacobian(dw,[wcx;wcy;wcz]));
partials_dwc=simplify(jacobian(dw,[dwcx;dwcy;dwcz]));
partials_k=simplify(jacobian(dw,[k1;k2]));


matlabFunction(partials_k,'File','STM_omega_k','Sparse',true);
matlabFunction(partials_w,'File','STM_omega_w','Sparse',true);
matlabFunction(partials_s,'File','STM_omega_s','Sparse',true);
% matlabFunction(partials_wc,'File','STM_omega_wc','Sparse',true);
% matlabFunction(partials_dwc,'File','STM_omega_dwc','Sparse',true);