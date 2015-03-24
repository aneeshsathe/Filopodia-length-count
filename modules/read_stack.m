function [ gray_img, raw_img] = read_stack( file_path )
%READ_STACK Summary of this function goes here
%   Detailed explanation goes here

img_info=imfinfo(file_path);


for count=1:numel(img_info)
    raw_img(:,:,count)=imread(file_path,'Index',count);
    
    gray_img(:,:,count)=mat2gray(raw_img(:,:,count));
end

end
