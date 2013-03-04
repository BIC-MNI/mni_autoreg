%
%  apply the non-linear transform in warps to the point p
%
function x = transf(p, warps);

t = inter_p( warps, p);

x = p+t;


