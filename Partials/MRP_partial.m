%% partials MRP
clear
clc
close all


syms s1 s2 s3 wx wy wz

s=[s1;s2;s3];
w=[wx;wy;wz];
skew_s=[0,-s3,s2;s3,0,-s1;-s2,s1,0];
C=eye(3)+8*(skew_s*skew_s)/(1+transpose(s)*s)^2-4*(1-transpose(s)*s)/(1+transpose(s)*s)^2*skew_s;
B=1/4*((1-transpose(s)*s)*eye(3)+2*skew_s+2*s*transpose(s));

ds=B*w;

partials=simplify(jacobian(ds,[wx;wy;wz;s1;s2;s3]));

matlabFunction(partials,'File','STM_MRP','Sparse',true);