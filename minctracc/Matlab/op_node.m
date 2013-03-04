function d = op_node(mod,md,mdd ,dat,d,dd, spacing, warps, pos,flag,op_type)

%   check to see if we are too close to the ends

if ((pos+spacing > max(md(:,2))) | (pos-spacing < min(md(:,2))))
   d = 0;
   return
end

%   get the value of the 1st deriv in the neighbourhood of the point of interest

for i=-4:4,
  p_i = pos + i*spacing/4;
  md_val = inter_1p(mod, p_i);
  if abs(md_val) > 0.2*max(md(:,1))
     break
  end;
end

%   if the value is not large enough, then ignore this node

if abs(md_val) < 0.2*max(md(:,1))
   d = 0;
   return;
end

vspace = spacing / (md(2,2)-md(2,1));

data = zeros(19,2);
model= zeros(11,2);

for i=-5:5,
  p_i = pos + i*spacing/4;
  model(i+6,1) = inter_2p(mod, p_i);
  model(i+6,2) = p_i;
end

%   get values of data 2nd derivative

for i=-9:9,
   pd = pos + i*spacing/4;
%				apply current transformation to get data point.
   p_i = transf(pd, warps);
   data(i+10,1) = inter_2p(dat, p_i);
   data(i+10,2) = pd;
end;

[d,tmpp] = match(model(:,:),data(:,:),spacing,op_type);

LR  = max(tmpp(:,2)) - min(tmpp(:,2));
DIFF= max(tmpp(:,1)) - min(tmpp(:,1));

if (flag>1)
 figure(5);
 plot(model(:,2),model(:,1),'x',data(:,2),data(:,1),'o'); hold on
 plot(model(:,2),model(:,1),data(:,2),data(:,1)); 
 if op_type==1
   plot([model(5,2) model(5,2)+d], [model(5,1) model(5,1)],'-w');
 else
   plot([model(5,2) model(5,2)+d], 0.8*[model(5,1) model(5,1)],'-g');
 end

 mt = max(abs(tmpp(:,1)));
 tmpp(:,1) = tmpp(:,1) / mt;
 tmpp(:,1) = max(abs(model(:,1))) * tmpp(:,1);
 plot(tmpp(:,2),tmpp(:,1),'r');
 grid on ; 
 hold off;
 title('model (x), data (o)');

 drawnow;
end;

