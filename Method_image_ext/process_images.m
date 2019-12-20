function [yn]=process_images(renderer,frame,method,cam_params)

        [J2_L,J2_R,disparityMap]=img_preproc(renderer,frame);
        
switch method
    case 'depth'
        [feats_HK, measures]=extract_features(J2_L,disparityMap,cam_params);
    otherwise
        [feats_HK, measures]=extract_and_match_features(J2_L,J2_R,cam_params);
end

        yn.m=size(measures,1);
        yn.z=double(reshape(measures',[],1)); % cast to double to use in functions
        yn.feats=feats_HK;



end