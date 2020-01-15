function  [features, filt_match1, filt_match2]=extract_and_match_features(I1gray,I1gray_R,cam_params)


% features detection:
blobs1H = detectHarrisFeatures(I1gray,'MinQuality',0.1,'FilterSize',5);
% blobs1H = detectSURFFeatures(I1gray);

[a,b]=sort(blobs1H.Metric,'desc');
sorted_blobs1_f=blobs1H(b);
% sorted_blobs1_f=selectStrongest(blobs1H,min(20,size(blobs1H,1)));

blobs2H = detectHarrisFeatures(I1gray_R,'MinQuality',0.1,'FilterSize',5);
% blobs2H = detectSURFFeatures(I1gray_R);

[a,b]=sort(blobs2H.Metric,'desc');
sorted_blobs2_f=blobs2H(b);
% sorted_blobs2_f=selectStrongest(blobs2H,min(20,size(blobs2H,1)));

% features extraction
[feats1HK, valid_points1_HK]=extractFeatures(I1gray,sorted_blobs1_f,'Method','FREAK','FeatureSize',128);
[feats2HK, valid_points2_HK]=extractFeatures(I1gray_R,sorted_blobs2_f,'Method','FREAK','FeatureSize',128);
% [feats1HK, valid_points1_HK]=extractFeatures(I1gray,sorted_blobs1_f,'Method','Block');
% [feats2HK, valid_points2_HK]=extractFeatures(I1gray_R,sorted_blobs2_f,'Method','Block');

% features matching
indexPairs = matchFeatures(feats1HK,feats2HK);
matchedPoints1_HK = valid_points1_HK(indexPairs(:,1),:);
matchedPoints2_HK = valid_points2_HK(indexPairs(:,2),:);

if  isa(feats1HK,'binaryFeatures')
matchedFeats1_HK=feats1HK.Features(indexPairs(:,1),:);
else
matchedFeats1_HK=feats1HK(indexPairs(:,1),:);
end

% additional filter as I know disparity only along horizontal lines
filter=abs(matchedPoints1_HK.Location(:,2)-matchedPoints2_HK.Location(:,2))<1;

filt_match1=matchedPoints1_HK(filter);
filt_match2=matchedPoints2_HK(filter);
filt_matchf1=matchedFeats1_HK(filter,:);

xq=filt_match1.Location(:,1);
yq=filt_match1.Location(:,2);

xqr=filt_match2.Location(:,1);

vq=xqr-xq;

measures=[xq-cam_params.hpix/2,-yq+cam_params.vpix/2,vq];

if  isa(feats1HK,'binaryFeatures')
    features=binaryFeatures(filt_matchf1);
else
    features=filt_matchf1;
end
    

% figure
% imshowpair(I1gray,I1gray_R,'montage');
% hold on
% plot(sorted_blobs1_f)
% plot(sorted_blobs2_f)
% 
% figure
% showMatchedFeatures(I1gray,I1gray_R,filt_match1,filt_match2);


end