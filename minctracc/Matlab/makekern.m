function kern = makekern(vsize,fwhm,size,type)

kern = zeros(size,1);

    
for k = -round(size/2):(round(size/2)-1),

  kindex = rem((k + size), size) +1;

  if (type=='g')
     kern(kindex) = normdist(vsize,vsize,fwhm,0.0,vsize*k);
  elseif (type=='r')
     kern(kindex) = rectdist(vsize,vsize,fwhm,0.0,vsize*k);
  end
end

if (type=='r')
   kern = kern / max(kern);
   c = ceil( fwhm/vsize);
   kern = kern ./ (fwhm/vsize +1);
   r = sum(kern);

   kern( round(c/2) + 1) =  kern( round(c/2) + 1) + (1 - r)/2;
   kern( size - round(c/2) +1 ) =  kern( size - round(c/2) +1 ) +  (1 - r)/2;
end
