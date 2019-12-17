function compatibility = compute_compatibility (prediction, observations)
%-------------------------------------------------------
%-------------------------------------------------------

% Compute individual distances
compatibility.d2 = zeros(observations.m, prediction.n);

for i = 1:observations.m
     indi=reshape(3*i-2+[0:2]',[],1);
     z = observations.z(indi);
     for j = 1:prediction.n
         indj=reshape(3*j-2+[0:2]',[],1);
         e = z - prediction.h(indj);
         C = prediction.HPH(indj,indj);
         compatibility.d2(i,j) = mahalanobis(e, C);
     end
end

% individual measure of compatibility is 0.99 percentile with 3 variables
% (u,v,d):

compatibility.ic = compatibility.d2 < chi2inv(0.99,3);
compatibility.candidates.features = find(sum(compatibility.ic, 1));
compatibility.candidates.observations = find(sum(compatibility.ic, 2))';

compatibility.AL = (sum (compatibility.ic, 2))';
compatibility.HS = prod(compatibility.AL + 1);

