function d = int_like(origdata, template)


d = zeros(size(template));

d(:,2) = template(:,2);

m = length(template);

for i=1:m,
   d(i,1) = inter_p(origdata,d(i,2));
end

