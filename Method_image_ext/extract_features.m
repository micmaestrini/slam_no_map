function [feats_HK, measures]=extract_features(J,disparityMap,cam_params)


% blobs = detectHarrisFeatures(J,'MinQuality',0.05);
blobs = detectSURFFeatures(J);

[~,b]=sort(blobs.Metric,'desc');
sorted_blobs=blobs(b);

strongest_blobs=sorted_blobs(1:min(200,size(sorted_blobs,1)));

% x=linspace(1,size(disparityMap,1),size(disparityMap,1))';
% y=linspace(1,size(disparityMap,2),size(disparityMap,2))';

% [X,Y]=meshgrid(x,y);

[feats_HK, valid_points_HK]=extractFeatures(J,strongest_blobs,'Method','KAZE');

xq=(valid_points_HK.Location(:,1));
yq=(valid_points_HK.Location(:,2));
% [Xq,Yq]=meshgrid(xq,yq);
% vq=diag(disparityMap(xq,yq));
vq = interp2(disparityMap,xq,yq);

feats_HK=feats_HK(~isnan(vq),:);
measures=[xq(~isnan(vq))-cam_params.hpix/2,-yq(~isnan(vq))+cam_params.vpix/2,-vq(~isnan(vq))];

% figure()
% imshow((disparityMap),[208;240])
% title('Disparity Map')
% colormap jet
% colorbar
% hold on
% plot(valid_points_HK);
% 
% figure()
% scatter3(xq(Vq>=206),yq(Vq>=206),Vq(Vq>=206))
% xlim([1,size(disparityMap,1)])
% ylim([1,size(disparityMap,2)])
% Vq = interp2(X,Y,V,Xq,Yq)

% figure()
% surf(X,Y,disparityMap','EdgeColor','None');
% zlim([206,240]);
% xlim([1,size(disparityMap,1)])
% ylim([1,size(disparityMap,2)])

% hold on
% 
% scatter3(xq,yq,239*ones(size(xq,1),1));
end