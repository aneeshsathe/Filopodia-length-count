function reshape_imshow( in_img, tool_or_show )
%RESHAPE_IMSHOW Summary of this function goes here
%   Detailed explanation goes here

%%
% in_img=BWao;
if strcmp(tool_or_show,'show')
    imshow(reshape(in_img, size(in_img,1),size(in_img,2)*size(in_img,3)))
elseif strcmp(tool_or_show,'tool')
    imtool(reshape(in_img, size(in_img,1),size(in_img,2)*size(in_img,3)))
end
end

