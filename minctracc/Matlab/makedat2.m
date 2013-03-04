function data = makedat2(len,i)

data = zeros(len,1);

if (i==1)
  data = set_rang(data,round(0.2*len),round(0.5*len), 0.75);   % model
  data = set_rang(data,round(0.515*len),round(0.533*len), 0.40);   % 
  data = set_rang(data,round(0.55*len),round(0.61*len), 0.85);   % 
else
  data = set_rang(data,round(0.2*len),round(0.48*len), 0.80);   % data
  data = set_rang(data,round(0.5*len),round(0.515*len), 0.60);   % 
  data = set_rang(data,round(0.54*len),round(0.58*len), 0.90);   %
end

return


if (i==1)
  data = set_rang(data,round(0.2*len),round(0.5*len), 0.75);   % 
  data = set_rang(data,round(0.53*len),round(0.56*len), 0.20);   % everything
  data = set_rang(data,round(0.6*len),round(0.67*len), 0.85);   % everything
else
  data = set_rang(data,round(0.2*len),round(0.45*len), 0.80);   % 
  data = set_rang(data,round(0.5*len),round(0.53*len), 0.50);   % everything
  data = set_rang(data,round(0.56*len),round(0.63*len), 0.90);   % everything
end