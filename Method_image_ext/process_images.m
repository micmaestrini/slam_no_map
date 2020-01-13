function [yn,frameLeftGray]=process_images(renderer,frame,method,cam_params)

        [frameLeftGray,frameRightGray,disparityMap]=img_preproc(renderer,frame,cam_params);
        
switch method
    case 'depth'
        [feats_HK, measures]=extract_features(frameLeftGray,disparityMap,cam_params);
    otherwise
        [feats_HK, measures]=extract_and_match_features(frameLeftGray,frameRightGray,cam_params);
end

        yn.m=size(measures,1);
        yn.z=double(reshape(measures',[],1)); % cast to double to use in functions
        yn.feats=feats_HK;



end