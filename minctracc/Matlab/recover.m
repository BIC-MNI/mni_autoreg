% function warps = recover(m,md,mdd ,dat,d,dd, spacing, deforms, flag)
% m   - model blurred intensity
% md  - model 1st derivative
% mdd - model 2nd derivative
% dat - data blurred intensity
% d   - data 1st derivative
% dd  - data 2nd derivative
% spacing - scale at which to recover warp
% deforms - existing deformation, already recovered.
% flag    - =1, do graph of match; =0, do nothing

function warps = recover(mod, md,mdd ,dat,d,dd, spacing, deforms, flag)

global offsets;
global def_field;
global field_count;

d = max(md(:,2)) - min(md(:,2));

n = floor(d/spacing);

%	warps will store the additional warp for this scale step
warps = zeros(n,2);

warps(1,2) = min(md(:,2)) + 0.5*(d-n*spacing);

for i=2:n,
   warps(i,2) = warps(1,2) + (i-1)*spacing;
end

%	total is the total global deformation to date.
total = int_like(deforms,warps);

%	current will be the current sum of total + this iteration's warp
current = total;


% 	warps1 will store the partial warp, calculated at sub-iteration
%	smooth will be a smoothed version of warps1

warps1 = warps; 
smooth = warps; 

fprintf('for %g spacing:\n',spacing);

iters = 20;
n = length(warps(:,1));

for j = 1:iters,

   if j<=3,
      op_type=1;
   else
      op_type=0;
   end

   for i=1:n,
      warps1(i,1) = op_node(mod,md,mdd ,dat,d,dd, spacing, current,warps(i,2),flag,op_type);
   end

   smooth(:,1) = warps(:,1) + 0.4*warps1(:,1);

%	now smooth out the warp:

   for i=2:n-1,
      warps(i,1)=0.2*smooth(i-1,1) + 0.6*smooth(i,1) + 0.2*smooth(i+1,1);
   end
   warps(1,1)=0.5*(smooth(1,1)+smooth(2,1));
   warps(n,1)=0.5*(smooth(n,1)+smooth(n-1,1));

%	build up the current global warp

   current(:,1) = total(:,1) + warps(:,1);

   field_part= int_like(current,offsets);
   field_part(:,1) = -field_part(:,1);
   def_field(:,field_count) = field_part(:,1); field_count=field_count+1;


   fprintf ('%g of %g: max %g   mean %g\n', j, iters, max(abs(warps1(:,1))), mean(abs(warps1(:,1))) );

if (flag>0),

  w = warpdata(mdd, field_part);

  figure(6)
  subplot(211);
  plot(offsets(:,2),offsets(:,1),'g--');hold on
  plot(warps(:,2),-current(:,1),'y'); 
  plot(warps(:,2),-current(:,1) + warps(:,1) - warps1(:,1),'y:'); 
  plot(warps1(:,2),-warps1(:,1),'w');
  plot(warps(:,2),-warps(:,1),'r');
  grid on; hold off
  drawnow;

  subplot(212);

  plot(dd(:,2),dd(:,1),'g--');hold on
  plot(w(:,2),w(:,1),'y');
  grid on; hold off
  drawnow;

end

end;


return



