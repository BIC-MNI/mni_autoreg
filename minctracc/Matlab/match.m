%
%  d = match(model,data,limit,op_type)
%  find the horizontal offset to best match model with data
%
%  if op_type = 1, do brute force, otherwise use secant
function [dist,p] = match(model,data,limit,op_type)

m = model;
d = data;

%m(:,1) = m(:,1) ./ max(abs(model(:,1)));
%d(:,1) = d(:,1) ./ max(abs(data(:,1)));

step = limit/10;



if op_type==1,
  [dist,p] = brute('differ',m,d,step,limit);
else
  [dist,p] = gold('differ',m,d,limit);
end


