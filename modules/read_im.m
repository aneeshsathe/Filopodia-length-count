function [filo_img, actin_img, raw_filo_img, raw_actin_img]= read_im(file_path)
%READ_IM Summary of this function goes here
%   Detailed explanation goes here

% file_path='Z:\IT Core Work\Naila_filopodia\140826_HeLaJW_GFPmyo10_CheLifeA  5 0.03tps with and wo ATR drug\140826_HeLaJW_GFPmyo10_CheLA  5 0.03tps.tif';

raw_filo_img=imread(file_path,'Index',1);
raw_actin_img=imread(file_path,'Index',2);

filo_img=mat2gray(raw_filo_img);
actin_img=mat2gray(raw_actin_img);


end

