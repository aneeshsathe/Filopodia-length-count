function out_pt  = an_nearfar_pt( bw,in_cent,min_or_max )
%AN_DISTOFBW Gives distance of all points in a BW img from a specific
%point. Gives closest or farthest depending on option
%   Detailed explanation goes here
if min_or_max=='min'
    [a,b]=find(bw==1);
% [d,p] = min((a-in_cent(2)).^2 + (b-in_cent(1)).^2);
[~,p] = min((a-in_cent(2)).^2 + (b-in_cent(1)).^2);
%  d = sqrt(d);
elseif min_or_max=='max'
    [a,b]=find(bw==1);
% [d,p] = max((a-in_cent(2)).^2 + (b-in_cent(1)).^2);
[~,p] = max((a-in_cent(2)).^2 + (b-in_cent(1)).^2);
%  d = sqrt(d);
end

out_pt=[b(p), a(p)];

end

