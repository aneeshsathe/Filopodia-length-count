function D1 = an_D1_D2_find( BW4, perim_C, perim_R )
%AN_D1_D2_FIND Summary of this function goes here
%   Detailed explanation goes here
D1 = bwdistgeodesic(BW4, perim_C, perim_R,'chessboard');

end

