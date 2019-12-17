function [stm]=STM(X,dt)
mu=3.986*1e14;
x=X(1);
y=X(2);
z=X(3);
vx=X(4);
vy=X(5);
vz=X(6);
r=X(7);
dr=X(8);
df=X(10);
wx=X(11);
wy=X(12);
wz=X(13);
s1=X(14);
s2=X(15);
s3=X(16);
k1=X(17);
k2=X(18);

wcx=X(19);
wcy=X(20);
wcz=X(21);
dwcx=X(22);
dwcy=X(23);
dwcz=X(24);





partials_nerm = STM_nerm(df,dr,mu,r,x,y,z);
partials_omega_w = STM_omega_w(k1,k2,s1,s2,s3,wcx,wcy,wcz,wx,wy,wz);
partials_omega_s = STM_omega_s(dwcx,dwcy,dwcz,k1,k2,s1,s2,s3,wcx,wcy,wcz,wx,wy,wz);
partials_omega_k = STM_omega_k(k1,k2,s1,s2,s3,wcx,wcy,wcz,wx,wy,wz);
partials_s = STM_MRP(s1,s2,s3,wx,wy,wz);

%% assemblaggio

P0=sparse(zeros(14));
P0(1:6,1:6)=partials_nerm;
P0(7:9,7:9)=partials_omega_w;
P0(7:9,10:12)=partials_omega_s;
P0(7:9,13:14)=partials_omega_k;
P0(10:12,7:12)=partials_s;


stm=P0*dt+sparse(eye(14));


end