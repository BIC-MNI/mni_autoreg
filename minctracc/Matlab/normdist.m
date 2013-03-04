function f = normal_dist(c,vsize,fwhm,mu,x)

sigma = fwhm/2.35482;

if (sigma==0) 
  if (x==mu)
     f = c;
  else
     f = 0;
  end
else 
  t1 = c / (sqrt(2*pi) * sigma);
  t2 = (x-mu)*(x-mu)/(2*sigma*sigma);
  t3 = exp(-t2);
  f = t1*t3;
end


