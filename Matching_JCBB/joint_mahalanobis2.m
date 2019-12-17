function [d2k, Hk, Ck, hk, zk] = joint_mahalanobis2 (prediction, observations, H)
%-------------------------------------------------------
%-------------------------------------------------------

% Compute joint distance for a hypothesis
[~,i, j] = find(H);

indi=reshape(3*i-2+[0:2]',[],1);
indj=reshape(3*j-2+[0:2]',[],1);

zk = observations.z(indi);
hk = prediction.h(indj);
% Rk = observations.R(indi,indi);
Ck = prediction.HPH(indj,indj);
Hk = prediction.H(indj,:);
d2k = mahalanobis (zk - hk, Ck);
