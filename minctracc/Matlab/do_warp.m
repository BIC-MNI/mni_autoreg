function [offsets,total,w24,w16,w8,w4] = do_warp(vsize, data_noise, disp_flag);

% ------------------- set up a spatial axis ------------------- 

global offsets;
global def_field;
global field_count;

x = [-64:63]' * vsize;

def_field=zeros(length(x),200);
field_count=1;

%-------------------  make a model data function ------------------- 

model = zeros(length(x),2);
model(:,1) = makedata2(length(x));
model(:,1) = model(:,1) / max(model(:,1));
model(:,2) = x;
model = blurdata(model,2*vsize);

% ------------------- make noisy data ------------------- 

rand('seed',sum(100*clock));
data = model;
data(:,1) = data(:,1) - data_noise/2 +  (data_noise*max(data(:,1)))*rand(size(model(:,1)));

% ------------------- warp the data so that it is different from the model

offsets = zeros(length(x),2);

offsets(:,2) = x;
offsets(:,1) = makeoffs(x,10);

warped = warpdata(data,offsets);

figure(1);
hold off;
clg;


plot(model(:,2),model(:,1), warped(:,2),warped(:,1));
ylabel('intensity');
xlabel('spatial position (mm)');
title('data and warped data');
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

mblur24 = blurdata(model, 24);
md24    = derivdat(mblur24);
mdd24   = derivdat(md24);

mblur16 = blurdata(model, 16);
md16    = derivdat(mblur16);
mdd16   = derivdat(md16);

mblur12  = blurdata(model, 12);
md12     = derivdat(mblur12);
mdd12    = derivdat(md12);

mblur8  = blurdata(model, 8);
md8     = derivdat(mblur8);
mdd8    = derivdat(md8);

mblur4  = blurdata(model, 4);
md4     = derivdat(mblur4);
mdd4    = derivdat(md4);


plot(blur24(:,2),blur24(:,1),d24(:,2),10*d24(:,1),dd24(:,2),40*dd24(:,1)); 
hold on; 
plot(mblur24(:,2),mblur24(:,1),'--',md24(:,2),10*md24(:,1),'--',mdd24(:,2),40*mdd24(:,1),'--'); 
grid on;
hold off;
title('model - -, data ---');  
ylabel('intensity');
xlabel('spatial position (mm)');

drawnow;

% ------------------- dewarp the data ------------------- 

first = zeros(size(offsets));
first(:,2) = offsets(:,2);

def_field(:,field_count) = first(:,1); field_count=field_count+1;

warps = recover(mblur24,md24,mdd24 ,blur24,d24,dd24, 12, first,disp_flag);
w24 = int_like(warps, offsets);
total = w24;
tot24 = total;
plot_res(1,24,offsets,tot24,model,warped,mblur24,blur24,mdd24,dd24);


warps = recover(mblur16,md16,mdd16 ,blur16,d16,dd16, 8, total,disp_flag);
w16 = int_like(warps, offsets);
total(:,1) = total(:,1) + w16(:,1);
tot16 = total;
plot_res(2,16,offsets,tot16,model,warped,mblur16,blur16,mdd16,dd16);

warps = recover(mblur12,md12,mdd12 ,blur12,d12,dd12, 6, total,disp_flag);
w12 = int_like(warps, offsets);
total(:,1) = total(:,1) + w12(:,1);
tot12 = total;
plot_res(3,12,offsets,tot12,model,warped,mblur12,blur12,mdd12,dd12);

warps = recover(mblur8,md8,mdd8 ,blur8,d8,dd8, 4, total,disp_flag);
w8 = int_like(warps, offsets);
total(:,1) = total(:,1) + w8(:,1);
tot8 = total;
plot_res(4,8,offsets,tot8,model,warped,mblur8,blur8,mdd8,dd8);


warps = recover(mblur4,md4,mdd4 ,blur4,d4,dd4, 2, total,disp_flag);
w4 = int_like(warps, offsets);
tot4 = total;
total(:,1) = total(:,1) + w4(:,1);

plot_res(5,4,offsets,tot4,model,warped,mblur4,blur4,mdd4,dd4);


% warps = recover(mblur4,md4,mdd4 ,blur4,d4,dd4, 2, total,disp_flag);
% w4_1 = int_like(warps, offsets);
% total(:,1) = total(:,1) + w4_1(:,1);

% warps = recover(mblur4,md4,mdd4 ,blur4,d4,dd4, 2, total,disp_flag);
% w4_2 = int_like(warps, offsets);
% total(:,1) = total(:,1) + w4_2(:,1);

% warps = recover(mblur4,md4,mdd4 ,blur4,d4,dd4, 2, total,disp_flag);
% w4_3 = int_like(warps, offsets);
% total(:,1) = total(:,1) + w4_3(:,1);

def_field(:,field_count)=-total(:,1); field_count = field_count+1;

for i = field_count:200,
  def_field(:,i)=offsets(:,1);
end;

save def_field

mi = min(min(def_field)); ma = max(max(def_field));

im = 64*(def_field - mi ) /(ma-mi);

figure(6); image(im);


return


for i=1:field_count
   s=sprintf('%g',i);
   plot(tot(:,2),def_field(:,i),tot(:,2),offs(:,1));
   text(-100,8,s); drawnow;
end

toff = offsets;
w_model=zeros(length(x),200);

for i=1:field_count
  fprintf('%g\n',i);
  toff(:,1) = def_field(:,i);
  w=warpdata(model, toff);

  w_model(:,i) = w(:,1);

end;



for i=1:field_count-4
   s=sprintf('%g',i);
   plot(model(:,2),w_model(:,i),warped(:,2),warped(:,1));
   text(0,0.1,s); drawnow;
end

figure(2); clg; subplot(211)
for i=1:field_count-4
   s=sprintf('%g',i);
   subplot(211); 
   plot(model(:,2),w_model(:,i),warped(:,2),warped(:,1));
   text(0,0.1,s); grid on; drawnow;
   subplot(212);	
   plot(tot(:,2),def_field(:,i),tot(:,2),offs(:,1));
   text(-100,8,s); grid on;drawnow;
end



figure(2); clg; subplot(211);
for i=1:field_count-4
   s=sprintf('%g',i);
   subplot(211); 

   plot(model(:,2),w_model(:,i),warped(:,2),warped(:,1));

   text(0,0.1,s); grid on; drawnow;
   subplot(212);	
   plot(tot(:,2),def_field(:,i),tot(:,2),offs(:,1));
   text(-100,8,s); grid on;drawnow;
end



w_mod = w_model;
for i=1:20,
   toff(:,1) = def_field(:,i);	
   w=warpdata(mdd24, toff);   
   w_model(:,i) = w(:,1);
end;
for i=21:40,
   toff(:,1) = def_field(:,i);	
   w=warpdata(mdd16, toff);   
   w_model(:,i) = w(:,1);
end;
for i=41:60,
   toff(:,1) = def_field(:,i);	
   w=warpdata(mdd12, toff);   
   w_model(:,i) = w(:,1);
end;
for i=61:80,
   toff(:,1) = def_field(:,i);	
   w=warpdata(mdd8, toff);   
   w_model(:,i) = w(:,1);
end;
for i=81:100,
   toff(:,1) = def_field(:,i);	
   w=warpdata(mdd4, toff);   
   w_model(:,i) = w(:,1);
end;



figure(2); clg; subplot(211);
for i=1:field_count-4
   s=sprintf('%g',i);
   subplot(211); 
   if i>0 & i<21,
      plot(model(:,2),w_model(:,i),dd24(:,2),dd24(:,1));
   elseif i>20 & i<41,
      plot(model(:,2),w_model(:,i),dd16(:,2),dd16(:,1));
   elseif i>40 & i<61,
      plot(model(:,2),w_model(:,i),dd12(:,2),dd12(:,1));
   elseif i>60 & i<81,
      plot(model(:,2),w_model(:,i),dd8(:,2),dd8(:,1));
   else i>80 & i<101,
      plot(model(:,2),w_model(:,i),dd4(:,2),dd4(:,1));
   end 
   text(0,0.1,s); grid on; drawnow;
   subplot(212);	
   plot(tot(:,2),def_field(:,i),tot(:,2),offs(:,1));
   text(-100,8,s); grid on;drawnow;
end
