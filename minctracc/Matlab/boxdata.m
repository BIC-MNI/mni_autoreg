function blur = boxdata(data,fwhm)

blur = zeros(size(data));
blur(:,2) = data(:,2);

vsize = blur(2,2) - blur(1,2);

b = makekern(vsize,fwhm,length(data(:,1)),'r');
B = fft(b);
D = fft(data(:,1));
blur(:,1) = real(ifft(B .* D));
