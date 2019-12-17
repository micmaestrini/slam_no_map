%% partials wc e dwc quat

clear
clc
close all

syms x y z vx vy vz r dr f df mu


ddf=-2*dr*df/r;
ddx=2*df*vy+ddf*y+df^2*x-mu*(r+x)/((r+x)^2+y^2+z^2)^1.5+mu/r^2;
ddy=-2*df*vx-ddf*x+df^2*y-mu*y/((r+x)^2+y^2+z^2)^1.5;
ddz=-mu*z/((r+x)^2+y^2+z^2)^1.5;
ddr=r*df^2-mu/r^2;

X=[x;y;z;r;f;vx;vy;vz;dr;df];
dX=[vx;vy;vz;dr;df;ddx;ddy;ddz;ddr;ddf];



d=sqrt(x^2+y^2+z^2);
b0=(x+d)/sqrt(2*d*(x+d));
b1=0;
b2=-z/sqrt(2*d*(x+d));
b3=y/sqrt(2*d*(x+d));


q=simplify([b0;b1;b2;b3]);
dq=simplify(jacobian(q,X)*dX);
ddq=simplify(jacobian(dq,X)*dX);

q_x=simplify(jacobian(q,[x;y;z;vx;vy;vz]));
dq_x=simplify(jacobian(dq,[x;y;z;vx;vy;vz]));
ddq_x=simplify(jacobian(ddq,[x;y;z;vx;vy;vz]));

% 
% matlabFunction(q,'File','b','Sparse',true);
% matlabFunction(dq,'File','db','Sparse',true);
% matlabFunction(ddq,'File','ddb','Sparse',true);
% matlabFunction(q_x,'File','STM_b_x','Sparse',true);
% matlabFunction(dq_x,'File','STM_db_x','Sparse',true);
% matlabFunction(ddq_x,'File','STM_ddb_x','Sparse',true);


