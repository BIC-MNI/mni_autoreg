function offsets = makeoffs(x,maxoff);

offsets = zeros(size(x));
d = (max(x) - min(x))/2;
m = length(offsets);

for i = 1:m,
   angle = pi*x(i)/d;    % angle will make one full cycle
   offsets(i) = maxoff*(0.4*sin(angle) - 0.2*sin(2*angle+pi/2) + 0.4*sin(5*angle+pi/2) + 0.1*sin(5*angle+pi/2));
end;
