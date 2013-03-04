function c = corrdata(d1,d2)

c = zeros(size(d1));
c(:,2) = d1(:,2);

m = length(d1); 

w1 = max(d1(:,2)) - min(d1(:,2));
w2 = max(d2(:,2)) - min(d2(:,2));
m1 = min(d1(:,2));
m2 = min(d2(:,2));
vsize = abs(d1(2,2)-d1(1,2));
p2 = d2;

start = -round(m/2);
done  = round(m/2);

for i=1:m,
i
	offset = m1 - m2  - w2/2.0 + (i-1)*vsize;
	p2(:,2) = d2(:,2) + offset;
	p1 = int_like(d1,p2);
	c(i,1) = corrdat(p1,p2);
end

