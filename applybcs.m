function [force,ifix,ubar] = applybcs(x,y,numnod,length,height,g)

force = zeros(1,numnod*2);
ifix = zeros(1,numnod*2);
ubar = zeros(1,numnod*2);

for i=1:numnod
   %left
   if (x(i) == 0)
      ifix(2*i-1) = 1.0;
      ubar(2*i-1) = g;
   end
   %right
   if (x(i) == length)
      ifix(2*i-1) = 1.0;
      ubar(2*i-1) = -g;
   end
   %top
   if (y(i) == height)
      ifix(2*i) = 1.0;
      ubar(2*i) = -g;
   end
   %bottom
   if (y(i) == 0)
      ifix(2*i) = 1.0;
      ubar(2*i) = g;
   end
end