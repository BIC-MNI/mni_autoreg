function [dist,p] = gold(diff_fun,m,d,limit)

p=zeros(3,2);
p(1,2) = (max(m(:,2))+min(m(:,2)))/2;
p(1,1) =  feval(diff_fun,m, d, 0);

p(2,2) = p(1,2) - 0.6*limit;
p(2,1) =  feval(diff_fun,m, d, -0.6*limit);

p(3,2) = p(1,2) + 0.6*limit;
p(3,1) =  feval(diff_fun,m, d, 0.6*limit);
j = 4;

dist = 0;

ax=-0.6*limit;
cx=0.6*limit;


R = 0.61803399;
C = 1.0-R;

x0=ax;
bx = 0;
x3=cx;

if abs(cx-bx) > abs(bx-ax),
	x1=bx;
	x2=bx+C*(cx-bx);
else
	x2=bx;
	x1=bx-C*(bx-ax);
end

f1=feval(diff_fun,m,d,x1);
f2=feval(diff_fun,m,d,x2);

while (abs(x3-x0) > limit/16),
	if f2 < f1,
		x0=x1; x1=x2; x2=R*x1+C*x3;
		f0=f1; f1=f2; f2=feval(diff_fun,m,d,x2);
		p(j,1)=f2;
		p(j,2) = p(1,2) + x2;
	else
		x3=x2; x2=x1; x1=R*x2+C*x0;
		f3=f2; f2=f1; f1=feval(diff_fun,m,d,x1);
		p(j,1)=f1;
		p(j,2) = p(1,2) + x1;
	end
end

if (f1 < f2),
	dist=x1;
else
	dist=x2;
end


