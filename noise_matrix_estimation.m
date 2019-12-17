clear
clc
close all

Vf=24;
Hf=36;

npix=100;
Hpix=npix*1.5;
Vpix=npix;
H_int=Hpix+1;
V_int=Vpix+1;

N=10000000;

x=Hf*(0.5-rand(N,1));
y=Vf*(0.5-rand(N,1));
v=[x,y];


edge_x_pix=(-Hpix/2:1:Hpix/2)*Hf/Hpix;
edge_y_pix=(-Vpix/2:1:Vpix/2)*Vf/Vpix;

xc_pix=(edge_x_pix(2:end)+edge_x_pix(1:end-1))/2;
yc_pix=(edge_y_pix(2:end)+edge_y_pix(1:end-1))/2;


xr=xc_pix(discretize(x,edge_x_pix))';
yr=yc_pix(discretize(y,edge_y_pix))';


vr=[xr,yr];

e=v-vr;
mean(e)
cov(e)

% 
% scatter(x,y)
% hold on
% scatter(xr,yr)
