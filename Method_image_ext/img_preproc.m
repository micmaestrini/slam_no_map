function [frameLeftGray,frameRightGray,disparityMap]=img_preproc(renderer,img_index,cam_params)


fx=cam_params.foc*cam_params.hpix/cam_params.Hf;
fy=cam_params.foc*cam_params.vpix/cam_params.Vf;
cam_1=cameraParameters('IntrinsicMatrix',[fx,0,0;0,fy,0;1024,972,1],'ImageSize',[1944,2048]);
cam_2=cameraParameters('IntrinsicMatrix',[fx,0,0;0,fy,0;1024,972,1],'ImageSize',[1944,2048]);
stereoParams = stereoParameters(cam_1,cam_2,eye(3),[1;0;0]);



currentDir=pwd;
folders_struct=split(currentDir,'\');

images_root=strjoin({folders_struct{1:end-1}},'\');




filename_L=strcat(images_root,'\images_',renderer,'\img_L');
filename_R=strcat(images_root,'\images_',renderer,'\img_R');
fft='.tif';

frameName_L=char(strcat(filename_L,num2str(img_index),fft));
frameName_R=char(strcat(filename_R,num2str(img_index),fft));


I_RGB_L=imread(frameName_L);
% I_L=rgb2gray(I_RGB_L);
% J_L=adapthisteq(I_L,'NumTiles',[64 64],'ClipLimit',0.1,'Distribution','rayleigh');
% J2_L=wiener2(J_L);
% J2_L=I_L;

% figure()
%     imshowpair(I_L,J2_L,'montage');
     
I_RGB_R=imread(frameName_R);
% I_R=rgb2gray(I_RGB_R);
% J_R=adapthisteq(I_R,'NumTiles',[64 64],'ClipLimit',0.1,'Distribution','rayleigh');
% J2_R=wiener2(J_R);
% J2_R=I_R;



[J1_rect,J2_rect] = rectifyStereoImages(I_RGB_L,I_RGB_R,stereoParams);


frameLeftGray = rgb2gray(J1_rect);
frameRightGray = rgb2gray(J2_rect);
    
    

% figure()
%     imshowpair(I_R,J2_R,'montage');
    %%
%     close all
disparityRange=[208, 240];
% disparityMap = disparityBM(wiener2(J2_L),wiener2(J2_R),'DisparityRange',disparityRange,'UniquenessThreshold',100,'BlockSize',5);
disparityMap = disparitySGM((frameLeftGray),(frameRightGray),'DisparityRange',disparityRange,'UniquenessThreshold',20);
disparityMap=medfilt2(disparityMap);

% figure()
% imshow((disparityMap),disparityRange)
% title('Disparity Map')
% colormap jet
% colorbar

disparityMap(disparityMap<208)=nan;
disparityMap(disparityMap>238.9)=nan;


% figure()
% imshow(J2_L);
% figure()
% imshow(J2_R);
% hold on
% plot(A);


end