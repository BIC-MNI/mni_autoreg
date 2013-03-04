function data = makedata(len)

data = zeros(len,1);

left = 0.20;
right= 0.80;

data = set_rang(data,round((left)*len),round((right)*len), 0.10);   % everything

data = set_rang(data,round((left)*len),round((left+0.02)*len), 0.83);   % scalp
data = set_rang(data,round((right-0.02)*len),round((right)*len), 0.80);   % scalp

data = set_rang(data,round((left+0.04)*len),round((left+0.07)*len), 0.60);   % marrow
data = set_rang(data,round((right-0.07)*len),round((right-0.04)*len), 0.65);  % marrow

data = set_rang(data,round((left+0.10)*len),round((right-0.10)*len), 0.90);   % brain

data = set_rang(data,round(0.40*len),round(0.48*len), 0.20);  % ventricle
data = set_rang(data,round(0.52*len),round(0.60*len), 0.20);  % ventricle

