%% partials NERM
clear
clc
close all

syms x y z vx vy vz r dr f df mu


ddf=-2*dr*df/r;

ddx=2*df*vy+ddf*y+df^2*x-mu*(r+x)/((r+x)^2+y^2+z^2)^1.5+mu/r^2;

ddy=-2*df*vx-ddf*x+df^2*y-mu*y/((r+x)^2+y^2+z^2)^1.5;

ddz=-mu*z/((r+x)^2+y^2+z^2)^1.5;

ddr=r*df^2-mu/r^2;

dX=[vx;vy;vz;ddx;ddy;ddz];


partials=simplify(jacobian(dX,[x;y;z;vx;vy;vz]));

matlabFunction(partials,'File','STM_nerm','Sparse',true);