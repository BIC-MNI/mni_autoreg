function warped = warpdata(data,offsets);

warped = data;

m = length(warped(:,1));

start = warped(1,2);
step  = warped(1,2) - warped(1,1);

for i = 1:m,
   p = warped(i,2);
   p = transf(p, offsets);
   warped(i,1) = inter_p(data, p);
end;
