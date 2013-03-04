function deriv = estderiv(data,vsize)

deriv = zeros(size(data));

for j=2:(length(data)-2),
    
   deriv(j) = -data(j-1)/(3*vsize)  -data(j)/(2*vsize) +data(j+1)/vsize -data(j+2)/(6*vsize)  ;
%   deriv(j) = data(j+1)  -data(j);

end

