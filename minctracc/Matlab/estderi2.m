function deriv2 = estderi2(data,vsize)

deriv2 = zeros(size(data));

for j=2:(length(data)-2),
    
   deriv2(j) = data(j-1)/vsize  -2*data(j)/vsize + data(j+1)/vsize;

end

