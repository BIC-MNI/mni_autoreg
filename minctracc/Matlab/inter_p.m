%
%  real space value interpolation
%
function v = inter_p(data,real_p)

start = data(1,2);
step  = data(2,2) - data(1,2);

pos = 1 + (real_p-start)/step;


v = inter(data, pos);


