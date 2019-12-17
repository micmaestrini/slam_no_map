function answer = jointly_compatible (prediction, observations, H)
%-------------------------------------------------------
% University of Zaragoza
% Centro Politecnico Superior
% Robotics and Real Time Group
% Authors:  J. Neira, J. Tardos
% Date   :  7-2004
%-------------------------------------------------------
%-------------------------------------------------------

d2 = joint_mahalanobis2 (prediction, observations, H);
dof = 3*length(find(H));
answer = d2 < chi2inv(0.99,dof);


