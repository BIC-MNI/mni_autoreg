function kern = make_kernel_FT(size,vsize)

kern = zeros(size,1);

f_sample_size = 1.0/(vsize*size);
factor = i * 2.0 * pi * f_sample_size;
    
for k = -round(size/2):(round(size/2)-1),
  kindex = rem((k + size), size) +1;
  kern(kindex) =  factor * k;
end
