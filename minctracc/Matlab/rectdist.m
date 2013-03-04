function f = rect_dist(c,vsize,fwhm,mu,x)

t = x-mu;

if (t >= -fwhm/2  & t <= fwhm/2  ) 
   f = c/fwhm;
else
   f = 0.0;
end


