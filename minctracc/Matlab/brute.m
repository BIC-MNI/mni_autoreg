function [dist,p] = brute(diff_fun,m,d,step,limit)

p=zeros(length([-limit:step:limit]),2);
p(1,2) = m(1,2);

xstep = (max(m(:,2))-min(m(:,2)))/21;

v_min = 1000000;
dist = 0;


j=1;

for i = -limit:step:limit,
   v = feval(diff_fun,m, d, i);
   p(j,2) = p(1,2) + (j-1)*xstep;
   p(j,1) = v; j=j+1;
	
   if (v<v_min)
      v_min = v;
      dist = i;
   end
end

