%
%  d = differ(model,data,offset)
%  return the difference between the overlapping area
%  of the two functions of model+offset wrt data.
%
function val = differ(model,data,offset)

val = 0;
minx = min(data(:,2));
maxx = max(data(:,2));

for i=1:length(model)
   t = model(i,2) + offset;  % new position
   if (t>=minx & t <=maxx)
      d = model(i,1) - inter_p(data, t);
%      val = val + abs(d);
      val = val + d*d;
   end
end