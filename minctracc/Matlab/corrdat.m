function c = corrdat(d1,d2)

s1 = sqrt ( sum( d1(:,1) .^2 ) );
s2 = sqrt ( sum( d2(:,1) .^2 ) );

if (~(s1 == 0.0) & ~(s2 == 0.0)) 
   c  = sum( d1(:,1) .* d2(:,1)) / (s1*s2);
elseif (s1 == 0.0 & ~(s2 == 0.0))
   c = 0.0;
elseif (~(s1 == 0.0) & s2 == 0.0)
   c = 0.0;
else
   c = 1.0;
end


