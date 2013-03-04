%   This script will produce graphs to show an example profile
%   blurred with different kernels, and differentiated.
%
%   The data has 5 per cent added noise.

% ------------------- set up a spatial axis ------------------- 

data_noise = 0.05;
vsize = 1;
x = [-128:127]' * vsize;

%-------------------  make a model data profile ------------------- 

model = zeros(length(x),2);
model(:,1) = makedata2(length(x));
model(:,1) = model(:,1) / max(model(:,1));
model(:,2) = x;
model = blurdata(model,2*vsize);

% ------------------- warp the data so that it is different from the model

offsets = zeros(length(x),2);
offsets(:,2) = x;
offsets(:,1) = makeoffs(x,8);

data = warpdata(model,offsets);

% ------------------- make noisy data ------------------- 

rand('seed',sum(100*clock));
warped = model;
warped(:,1) = warped(:,1) - data_noise/2 +  (data_noise*max(warped(:,1)))*rand(size(model(:,1)));


figure(1);
hold off;
clg;
subplot(4,1,1);

plot(warped(:,2),warped(:,1));
ylabel('intensity');
xlabel('spatial position (mm)');
title('original data');
grid on;
drawnow;



% ------------------- blur the data with a Gaussian kernel and get the derivatives

blur24 = blurdata(warped, 24);
d24    = derivdat(blur24);
dd24   = derivdat(d24);

blur16 = blurdata(warped, 16);
d16    = derivdat(blur16);
dd16   = derivdat(d16);

blur12  = blurdata(warped, 12);
d12     = derivdat(blur12);
dd12    = derivdat(d12);

blur8  = blurdata(warped, 8);
d8     = derivdat(blur8);
dd8    = derivdat(d8);

blur4  = blurdata(warped, 4);
d4     = derivdat(blur4);
dd4    = derivdat(d4,vsize);



subplot(4,1,2);


plot(warped(:,2),blur4(:,1), 'g:'); hold
plot(warped(:,2),blur8(:,1), 'r-.')
plot(warped(:,2),blur16(:,1),'m--')
plot(warped(:,2),blur24(:,1),'w-')

text(100, 0.71, '... 4mm');
text(100, 0.64, '-.- 8mm');
text(100, 0.57, '- - 16mm');
text(100, 0.50, '--- 24mm');


ylabel('intensity');
xlabel('spatial position (mm)');
title('blurred data at different scales');
grid on;
drawnow;

subplot(4,1,3);


plot(warped(:,2),d4(:,1), 'g:'); hold;
plot(warped(:,2),d8(:,1), 'r-.')
plot(warped(:,2),d16(:,1),'m--')
plot(warped(:,2),d24(:,1),'w-')

text(100, 0.17, '... 4mm');
text(100, 0.13, '-.- 8mm');
text(100, 0.09, '- - 16mm');
text(100, 0.05, '--- 24mm');


ylabel('intensity');
xlabel('spatial position (mm)');
title('first derivative data at different scales');
grid on;
drawnow;

subplot(4,1,4);

plot(warped(:,2),dd4(:,1) * 4 , 'g:'); hold;
plot(warped(:,2),dd8(:,1) * 8, 'r-.')
plot(warped(:,2),dd16(:,1)*16 ,'m--')
plot(warped(:,2),dd24(:,1)*24 ,'w-')

text(100, 0.17, '... 4mm');
text(100, 0.13, '-.- 8mm');
text(100, 0.09, '- - 16mm');
text(100, 0.05, '--- 24mm');


ylabel('intensity * scale');
xlabel('spatial position (mm)');
title('second derivative data at different scales');
grid on;
drawnow;

figure(2);
clg
hold off

x4 = [-126:3:127];
template = zeros(length(x4),2);
template(:,2) = x4';
b4 = boxdata(blur8, 4);
b8_4 = int_like(b4, template);
t = int_like(b8_4, blur8);

plot(blur8(:,2), blur8(:,1), 'y'); hold;
plot(t(:,2), t(:,1), 'r'); 

plot(t(:,2), blur8(:,1)-t(:,1), 'g'); 

grid on




