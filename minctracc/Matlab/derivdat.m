function deriv = derivdat(data,vsize)

deriv = zeros(size(data));
deriv(:,2) = data(:,2);

vsize = deriv(2,2) - data(1,2);

b = makeft(length(data),vsize);
D = fft(data(:,1));
deriv(:,1)= real(ifft(b .* D));
