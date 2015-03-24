function out_img = an_gray(I)
%AN_GRAY returns graythresholded img
%   Detailed explanation goes here

out_img=im2bw(I,graythresh(I));
end

