%% =====READ IMAGE=======
clear
clc
file_path=[pwd '\140826_HeLaJW_GFPmyo10_CheLifeA  5 0.03tps with and wo ATR drug\140826_HeLaJW_GFPmyo10_CheLA  5 0.03tps.tif'];


filo_img=mat2gray(imread(file_path,'Index',1));
actin_img=mat2gray(imread(file_path,'Index',2));
% filo_img=imread(file_path,'Index',1);
% actin_img=imread(file_path,'Index',2);


%% texture analysis

Efilo=mat2gray(entropyfilt(filo_img));
Eactin=mat2gray(entropyfilt(actin_img));

bw_actin=im2bw(Eactin,graythresh(Eactin));
BWao_actin=bwareaopen(bw_actin,2000);

nhood = true(9);
closeBWao_actin = imclose(BWao_actin,nhood);
imshow([Eactin bw_actin BWao_actin closeBWao_actin])
%%

bw_filo=im2bw(Efilo,graythresh(Efilo));
imtool([filo_img Efilo bw_filo bw_actin bw_filo&bw_actin],[])
% imtool([filo_img Efilo actin_img Eactin],[])

%%


filo_thresh=filo_img>graythresh(filo_img);
act_thresh=actin_img> graythresh(actin_img);

imshow([filo_thresh act_thresh])
% 
% or_img=filo_thresh|act_thresh;
% and_img=filo_thresh&act_thresh;
% 
% or_and_or_img=xor(or_img,and_img);
% 
% imshow([filo_img filo_thresh actin_img act_thresh;
%     or_img and_img or_and_or_img zeros(512,512)])

%%
% se = strel('disk',15);
% afterOpening = imopen(actin_img,se);
% imtool(afterOpening)

%%

H = fspecial('prewitt');
filt_img = imfilter(actin_img,H,'replicate');
imtool([filt_img actin_img]);






