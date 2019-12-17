function [Pn, Kn, Hn]=propagate_aftermatch(X0,Xn,S0,P0,R0,Q0,match,foc,b,mu)
%% part where subset of points are identified:
% match array contains 2 columns, in the first one there are matched
% propagated points index. In the second column there are matched
% propagated measures.
%% propagate covariance:

vis=match(:,1);
n=size(match,1);
Sn=S0(match(:,1),:);

[Pn]=prop_covariance(X0,P0,Q0,vis,mu);
H_=sparse(5*n,15+3*n);
q0=Xn(14:17);
x0=Xn(1);
y0=Xn(2);
z0=Xn(3);

for i=1:n
    h_xi=H(b,foc,q0(1),q0(2),q0(3),q0(4),Sn(i,1),Sn(i,2),Sn(i,3),x0,y0,z0);
    h_pi=H_pi(b,foc,q0(1),q0(2),q0(3),q0(4),Sn(i,1),Sn(i,2),Sn(i,3),x0,y0,z0);
    H_((i-1)*5+1:i*5,1:6)=h_xi(:,1:6);
    H_((i-1)*5+1:i*5,10:13)=h_xi(:,7:end);
    H_((i-1)*5+1:i*5,15+3*i-2:15+3*i)=h_pi;   
end

Zn=H_*Pn*H_'+R0;

Kn=Pn*H_'/Zn;

Hn=H_;

end