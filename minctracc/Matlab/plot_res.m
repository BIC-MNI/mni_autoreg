function p=plot_res(fig,res,offsets,total,model,data,bmodel,bdata,mdd,dd);
figure(fig);

tot=total;
tot(:,1)=-1*total(:,1);

subplot(2,2,1)
plot(tot(:,2),tot(:,1)); hold on;
plot(offsets(:,2),offsets(:,1),'g--');
d=int_like(tot,offsets);
plot(offsets(:,2),abs(offsets(:,1)-d(:,1)),'w-.');
grid on;
hold off;
s= sprintf('%gmm:  true - -, recov ---, diff -.-',res); 
title(s);  
ylabel('warp offset (mm)');
xlabel('spatial position (mm)');

subplot(2,2,2)
plot(model(:,2),model(:,1),'y--'); hold on;
plot(data(:,2),data(:,1),'g-.');
d = warpdata(model,tot);
plot(d(:,2),d(:,1),'w');
grid on;
hold off;
title('model - -, orig -.-., recov ---');
ylabel('pixel intensity');
xlabel('spatial position (mm)');
drawnow;

subplot(2,2,3)
plot(bmodel(:,2),bmodel(:,1),'y--'); hold on;
plot(bdata(:,2),bdata(:,1),'g-.');
d = warpdata(bmodel,tot);
plot(d(:,2),d(:,1),'w');
grid on;
hold off;
title('model - -, orig -.-., recov ---');
ylabel('blurred intensity');
xlabel('spatial position (mm)');
drawnow;

subplot(2,2,4)
plot(mdd(:,2),mdd(:,1),'y--'); hold on;
plot(dd(:,2),dd(:,1),'g-.');
d = warpdata(mdd,tot);
plot(d(:,2),d(:,1),'w');
grid on;
hold off;
title('m - -, d -.-., recov ---');
ylabel('2nd deriv intensity');
xlabel('spatial position (mm)');
drawnow;

