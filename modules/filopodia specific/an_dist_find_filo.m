% function [ path_length ] = an_dist_find_filo( D2,D1)
function [ path_length ] = an_dist_find_filo( BW4,perim_C, perim_R,branch_C, branch_R )
%AN_DIST_FIND_FILO Summary of this function goes here
%   Detailed explanation goes here

D1 = bwdistgeodesic(BW4, perim_C, perim_R,'chessboard');
D2=bwdistgeodesic(BW4, branch_C, branch_R,'chessboard');

%     imshow(cell_perim)
%     hold on
%     scatter(branch_C(1200), branch_R(1200))
D = D1 + D2;
D = round(D * 8) / 8;

D(isnan(D)) = inf;
skeleton_path = imregionalmin(D);

pro_path_length = D(skeleton_path);
% disp('hello')
path_length= pro_path_length(1);
%         chess_dist=path_length;

end

