function [bw_in_img, img_props] = thresh_and_centroid( in_img )
%THRESH_AND_CENTROID Thresholds input img with graythresh and returns
%centroid of largest object

bw_in_img=im2bw(mat2gray(in_img),graythresh(in_img));


bw_in_img=keep_max_area_obj( bw_in_img);
bw_in_img=imfill(bw_in_img,'holes');

img_props=regionprops(bw_in_img);


% out_arg=floor(img_props.Centroid) ;

end

