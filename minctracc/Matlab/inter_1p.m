%
%  real space 2nd derivative value interpolation
%
function v = inter_1p(data,real_p)

start = data(1,2);
step  = data(2,2) - data(1,2);

%				establish voxel position
%				note: voxels start at 1
pos = 1 + (real_p-start)/step;
m = length(data);

p = floor(pos);
u = pos-p;

v=0;

if (pos>=2 &  pos<(m-1))
   v0 = data(p-1);
   v1 = data(p);
   v2 = data(p+1);
   v3 = data(p+2);

   v = cubic_p(v0,v1,v2,v3,u) / step;

elseif (pos>=1 &  pos<2)
   v0 = data(p);
   v1 = data(p+1);
   v2 = data(p+2);
   v3 = data(p+3);
   u=u+1;

   v = cubic_p(v0,v1,v2,v3,u) / step;

elseif (pos>=(m-1) &  pos<m)
   v0 = data(p-2);
   v1 = data(p-1);
   v2 = data(p);
   v3 = data(p+1);
   u=u-1;

   v = cubic_p(v0,v1,v2,v3,u) / step;

end
