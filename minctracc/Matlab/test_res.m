function  test_res(data, blur_size, voxel_size)

% this function will take the input data and blur it with a Gaussian kernel
% with a fwhm = blur_size.  The first and second derivitives will be
% calculated (in the Fourier domain) from the blurred data.



% ------------------- blur the data with a Gaussian kernel and get the derivatives

vsize = data(2,2) - data(1,2);

b = blurdata(data, blur_size);
d = derivdat(b);
dd= derivdat(d);


subplot(4,2,1);
plot(data(:,2),data(:,1),'g-')
ylabel('intensity');
xlabel('spatial position (mm)');
s = sprintf('original data (voxel size = %d mm)',vsize);
title(s);
grid on;

subplot(4,2,3);
plot(b(:,2),b(:,1),'g-');
ylabel('intensity');
xlabel('spatial position (mm)');
s = sprintf('blurred w/gaussian (fwhm=%4.1f, voxel=%4.1f)',blur_size,voxel_size);
title(s);
grid on;

subplot(4,2,5);
plot(d(:,2),d(:,1),'g-');
ylabel('intensity');
xlabel('spatial position (mm)');
s = sprintf('first derivative ');
title(s);
grid on;

subplot(4,2,7);
plot(dd(:,2),dd(:,1),'g-');
ylabel('intensity');
xlabel('spatial position (mm)');
s = sprintf('second derivative ');
title(s);
grid on;


% decimate the blurred data with the required voxel size


start = min(data(:,2)) - vsize/2 + voxel_size/2;
last  = max(data(:,2)) + vsize/2 - voxel_size/2;

x = [start:voxel_size:last];
template = zeros(length(x),2);
template(:,2) = x';
btmp = boxdata(b, voxel_size);
sb = int_like(btmp, template);
tb = int_like(sb, b);

subplot(4,2,3); hold on;
plot(tb(:,2), tb(:,1), 'y--');
plot(tb(:,2), b(:,1)-tb(:,1), 'r-.'); 

% using the decimated data, estimate the 1st derivative, using
% cubic interpolation

td = d;
td(:,1) = zeros(size(td(:,1)));
for i=1:length(td)
	td(i,1) = inter_1p( sb, td(i,2) );
end

%td = int_like(sd, d);
subplot(4,2,5); hold on;
plot(td(:,2), td(:,1), 'y--');
plot(td(:,2), d(:,1)-td(:,1), 'r-.'); 

% using the decimated data, estimate the 2nd derivative, using
% cubic interpolation

tdd = dd;
tdd(:,1) = zeros(size(tdd(:,1)));
for i=1:length(tdd)
	tdd(i,1) = inter_2p( sb, tdd(i,2) );
end

%tdd = int_like(sdd, dd);
subplot(4,2,7); hold on;
plot(tdd(:,2), tdd(:,1), 'y--');
plot(tdd(:,2), dd(:,1)-tdd(:,1), 'r-.'); 


