function [frameLeftGray,filt_match1,filt_match2, feats_HK]=process_images(renderer,frame,method,cam_params)

        [frameLeftGray,frameRightGray,disparityMap]=img_preproc(renderer,frame,cam_params);
        
switch method
    case 'depth'
        [feats_HK, filt_match1, filt_match2]=extract_features(frameLeftGray,disparityMap,cam_params);
    otherwise
        [feats_HK, filt_match1, filt_match2]=extract_and_match_features(frameLeftGray,frameRightGray,cam_params);
end

end