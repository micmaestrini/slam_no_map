function [J2_L,J2_R,disparityMap]=img_preproc(renderer,img_index)

currentDir=pwd;
folders_struct=split(currentDir,'\');

images_root=strjoin({folders_struct{1:end-1}},'\');



filename_L=strcat(images_root,'\images_',renderer,'\img_L');
filename_R=strcat(images_root,'\images_',renderer,'\img_R');
fft='.tif';

frameName_L=char(strcat(filename_L,num2str(img_index),fft));
frameName_R=char(strcat(filename_R,num2str(img_index),fft));


I_RGB_L=imread(frameName_L);
I_L=rgb2gray(I_RGB_L);
J_L=adapthisteq(I_L,'NumTiles',[64 64],'ClipLimit',0.1,'Distribution','rayleigh');
J2_L=wiener2(J_L);
    

% figure()
    imshowpair(I_L,J2_L,'montage');
     
I_RGB_R=imread(frameName_R);
I_R=rgb2gray(I_RGB_R);
J_R=adapthisteq(I_R,'NumTiles',[64 64],'ClipLimit',0.1,'Distribution','rayleigh');
J2_R=wiener2(J_R);

% figure()
%     imshowpair(I_R,J2_R,'montage');
    %%
%     close all
disparityRange=[208, 240];
% disparityMap = disparityBM(wiener2(J2_L),wiener2(J2_R),'DisparityRange',disparityRange,'UniquenessThreshold',100,'BlockSize',5);
disparityMap = disparityBM((J2_L),(J2_R),'DisparityRange',disparityRange,'UniquenessThreshold',30,'BlockSize',35);



disparityMap=medfilt2(disparityMap);
% figure()
% imshow((disparityMap),disparityRange)
% title('Disparity Map')
% colormap jet
% colorbar
% hold on
% plot(A);


end