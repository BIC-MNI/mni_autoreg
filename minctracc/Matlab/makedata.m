function data = makedata(len)

data = zeros(len,1);



data = set_rang(data,round(0.08*len),round(0.92*len), 0.10);   % everything

data = set_rang(data,round(0.08*len),round(0.10*len), 0.83);   % scalp
data = set_rang(data,round(0.90*len),round(0.92*len), 0.80);   % scalp

data = set_rang(data,round(0.12*len),round(0.15*len), 0.60);   % marrow
data = set_rang(data,round(0.85*len),round(0.88*len), 0.65);  % marrow

data = set_rang(data,round(0.18*len),round(0.82*len), 0.90);   % brain

data = set_rang(data,round(0.40*len),round(0.48*len), 0.20);  % ventricle
data = set_rang(data,round(0.52*len),round(0.60*len), 0.20);  % ventricle

