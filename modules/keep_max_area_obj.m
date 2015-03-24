function [out_BW, BW] = keep_max_area_obj(BW)
%KEEP_MAX_AREA_OBJ returns BW img with only the object of max size from
%input BW img
%%
%  BW=an_gray(MCimg);


CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
to_rem=cat(1,CC.PixelIdxList{[1:idx-1,idx+1:size(CC.PixelIdxList,2)]});
BW(to_rem) = 0;

BW=imfill(BW,'holes');

se = strel('diamond', 3);
erBW = imopen(BW,se);

CC = bwconncomp(erBW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
to_rem=cat(1,CC.PixelIdxList{[1:idx-1,idx+1:size(CC.PixelIdxList,2)]});
erBW(to_rem) = 0;


out_BW=erBW;
 

% imshow([BW,erBW,an_gray(MCimg),(erBW+an_gray(MCimg)./2)])

end

