function makefig()

x = [-74:53] * 0.5;

d1 = zeros( length(x), 2);

d1(:,2) = x';
d2 = d1;

d1(:,1) = makedat2(length(x),1);
d2(:,1) = makedat2(length(x),2);


b1 = blurdata(d1,3);
b2 = blurdata(d2,3);

d1_d  = derivdat(b1);
d1_dd = derivdat(d1_d);
d2_d  = derivdat(b2);
d2_dd = derivdat(d2_d);


part = [50:100];

clg;

subplot(3,1,1);
plot(b1(part,2),b1(part,1),'y',b2(part,2),b2(part,1),'w--');   grid on;
title('model ---,  data - - -');
ylabel('image intensity');
xlabel('offset (mm)');

subplot(3,1,2);
plot(d1_d(part,2),3*abs(d1_d(part,1)),'y',d2_d(part,2),3*abs(d2_d(part,1)),'w--');grid on;
title('model ---,  data - - -');
ylabel('1st deriv ');
xlabel('offset (mm)');

subplot(3,1,3);
plot(d1_dd(part,2),3*d1_dd(part,1),'y',d2_dd(part,2),3*d2_dd(part,1),'w--');grid on;
title('model ---,  data - - -');
ylabel('2nd deriv ');
xlabel('offset (mm)');



