function out_BW = area_thresh_obj( BW, area_thresh )
%AREA_THRESH_OBJ Summary of this function goes here
%   Detailed explanation goes here
%%
% BW=remove_perim;
% area_thresh=3;
BW=logical(BW);
CC = regionprops(BW,'Area','PixelIdxList');

CC([CC(:).Area]<area_thresh)=[];

out_BW=zeros(size(BW));

for count=1:numel(CC)
out_BW([CC(count).PixelIdxList])=1;

end
% imshow(out_BW)
% numPixels = cellfun(@numel,CC.PixelIdxList);
% [~,idx] = max(numPixels);
% to_rem=cat(1,CC.PixelIdxList{[1:idx-1,idx+1:size(CC.PixelIdxList,2)]});
% BW(to_rem) = 0;
% 
% BW=imfill(BW,'holes');

end
