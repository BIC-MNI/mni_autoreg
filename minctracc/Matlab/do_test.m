function [ee,dd] = do_test(vsize, width, fwhm, data_noise, noise);

% make a unit height rect function

%data = fftshift( makekern(vsize, width, 128, 'r'));
data = makedata(128);
rand('seed',sum(100*clock))
data = data -data_noise/2 +  data_noise*max(data)*rand(size(data));
data = data / max(data);

% set up a spatial axis

x = [-64:63] * vsize;

% ----------------------------------

% blur the data with a Gaussian kernel 

blur = blurdata(data, fwhm, vsize);

% add noise

rand('seed',sum(100*clock))
noisy = blur -noise/2 +  noise*max(blur)*rand(size(blur));

% get the first derivative

d = derivdat(noisy,vsize);

% get the second derivative

dd = derivdat(d,vsize);

% plot the results
a = [1:128];
%a = [40:90];

plot(x(a),blur(a),x(a),data(a),x(a),d(a)*3,x(a),dd(a)*20, x(a), noisy(a));
grid on
xlabel('spatial axis, mm');
ylabel('arbitrary intensity');

fprintf(1,'hit a key\n');

% pause

% e = estderiv(noisy,vsize);
% ee= estderi2(noisy,vsize);
% plot(x(a),noisy(a),x(a),d(a),x(a),e(a),'--',x(a),5*dd(a),x(a),10*ee(a),'--');
% title ('dotted lines are cubic spline interpolation');
% xlabel('spatial axis, mm');
% ylabel('arbitrary intensity');

% grid on


