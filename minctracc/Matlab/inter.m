%
%  voxel space value interpolation, using
%  cubic where possible, linear near ends 
%  and nearest at ends.
%
function v = inter(data,pos)

m = length(data);

if (pos>=2 &  pos<(m-1))
   p = floor(pos);
   u = pos-p;

   v0 = data(p-1);
   v1 = data(p);
   v2 = data(p+1);
   v3 = data(p+2);

   v = ( (v1) + (u) * ( 0.5 * ((v2)-(v0)) + (u) * ( (v0) - 2.5 * (v1) + 2.0 * (v2) - 0.5 * (v3) + (u) * ( -0.5 * (v0) + 1.5 * (v1) - 1.5 * (v2) + 0.5 * (v3)  ) ) ) );

else
   if (pos<0.5 | pos>=(m+0.5))
      v = 0;
   else
      if (pos<1 | pos>m)
         v = 0.5 * data(round(pos));
      else
         if (pos<2)
            f = pos-1;
            v = (1.0-f)*data(1) + f*data(2);
         end
         if (pos>=m-1)
            f = pos-(m-1);
            v = (1.0-f)*data(m-1) + f*data(m);
         end
      end
   end
end
