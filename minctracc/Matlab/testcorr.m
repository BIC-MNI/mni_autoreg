x = [-128:127]' .* 0.5; 
y = zeros(length(x),2);
y(:,2) = x;
y(:,1) = [zeros(1,length(x)/8) 0.7*ones(1,length(x)/8) ones(1,length(x)/2) zeros(1,length(x)/4)]';


b = blurdata(y,width);
d = derivdat(b);

figure(1); clg;
subplot(4,1,1);
v = [1:128]';
plot(b(v,2),b(v,1)./20, d(v,2),d(v,1));


e = d(v,:);

f = [4:80]; 
p = zeros(size(f),2);
p = e(f,:);
c = corrdata(e,p);
subplot(4,1,2); 
plot(p(:,2),p(:,1).*10,'c--'); hold on;
plot(c(:,2),c(:,1),'y'); hold off

f = [23:42]; 
p = zeros(size(f),2);
p = e(f,:);
fc = corrdata(e,p); 
subplot(4,1,3); 
plot(p(:,2),p(:,1).*10,'c--'); hold on;
plot(fc(:,2),fc(:,1),'y'); hold off



f = [56:75]; 
p = zeros(size(f),2);
p = e(f,:);
fc = corrdata(e,p); 
subplot(4,1,4); 
plot(p(:,2),p(:,1).*10,'c--'); hold on;
plot(fc(:,2),fc(:,1),'y'); hold off


