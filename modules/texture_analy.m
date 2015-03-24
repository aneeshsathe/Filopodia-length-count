function [BWao, Eim]= texture_analy( I )
%TEXTURE_ANALY Gives texture based thresholding
%%
E = entropyfilt(I);
% imtool(E,[])
Eim = mat2gray(E);
%  imtool(Eim,[])
%%
BWao=zeros(size(Eim));
for count=1:size(Eim,3)
    BW1 = im2bw(Eim(:,:,count), .7);
    
    BWao(:,:,count) = bwareaopen(BW1,200);
    
end
% imshow(BWao);
% imtool([I Eim BW1 BWao],[])
% nhood = true(9);
% closeBWao = imclose(BWao,nhood);
% imshow(closeBWao)


end

