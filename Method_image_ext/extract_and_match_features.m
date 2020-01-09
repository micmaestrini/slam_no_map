function  [features, measures]=extract_and_match_features(I1gray,I1gray_R,cam_params)


% features detection:
blobs1H = detectHarrisFeatures(I1gray);
% blobs1H = detectSURFFeatures(I1gray);

[a,b]=sort(blobs1H.Metric,'desc');
sorted_blobs1=blobs1H(b);
sorted_blobs1=selectStrongest(blobs1H,min(100,size(blobs1H,1)));

blobs2H = detectHarrisFeatures(I1gray_R);
% blobs2H = detectSURFFeatures(I1gray_R);

[a,b]=sort(blobs2H.Metric,'desc');
sorted_blobs2=blobs2H(b);
sorted_blobs2=selectStrongest(blobs2H,min(100,size(blobs2H,1)));

% features extraction
[feats1HK, valid_points1_HK]=extractFeatures(I1gray,sorted_blobs1,'Method','FREAK','FeatureSize',128);
[feats2HK, valid_points2_HK]=extractFeatures(I1gray_R,sorted_blobs2,'Method','FREAK','FeatureSize',128);

% features matching
indexPairs = matchFeatures(feats1HK,feats2HK);
matchedPoints1_HK = valid_points1_HK(indexPairs(:,1),:);
matchedPoints2_HK = valid_points2_HK(indexPairs(:,2),:);

matchedFeats1_HK=feats1HK.Features(indexPairs(:,1),:);

% additional filter as I know disparity only along horizontal lines
filter=abs(matchedPoints1_HK.Location(:,2)-matchedPoints2_HK.Location(:,2))<1;

filt_match1=matchedPoints1_HK(filter);
filt_match2=matchedPoints2_HK(filter);
filt_matchf1=matchedFeats1_HK(filter,:);

xq=filt_match1.Location(:,1);
xqr=filt_match2.Location(:,1);
yq=filt_match2.Location(:,2);
vq=xqr-xq;

measures=[xq-cam_params.hpix/2,-yq+cam_params.vpix/2,vq];
features=binaryFeatures(filt_matchf1);

figure
showMatchedFeatures(I1gray,I1gray_R,filt_match1,filt_match2);


end