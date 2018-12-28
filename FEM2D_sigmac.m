function [xc,yc] = FEM2D_sigmac(lambda,mu,sigmac,length,height,Nx,Ny)
%INPUT DATA
[x,y,node,numele,numnod] = mesh2d(length,height,Nx,Ny);

%MATERIAL CONSTANTS
%already passed in
%a field, possibly random, of lambda
%a field, possibly random, of mu

%FORCE AND DISPLACEMENT BC'S
g = 1.0;
[force,ifix,ubar] = applybcs_sigmac(x,y,numnod,length,height,g);

%ASSEMBLY OF STIFFNESS
ndof = 2; %degrees of freedom per node
gauss = [-3^(-0.5), 3^(-0.5)];
numeqns = numnod*ndof;
bigk = zeros(numeqns);

%
% loop over elements
%
% nlink is # of nodes per element
nlink = 4;
for e = 1:numele
   [ke] = elemstiff(node,x,y,gauss,lambda,mu,e);
   %
   % assemble ke into bigk
   %
   n1 = ndof-1;
   for i=1:nlink
      for j=1:nlink
         rbk = ndof*(node(i,e)-1) + 1;
         cbk = ndof*(node(j,e)-1) + 1;
         re = ndof*(i-1)+1;
         ce = ndof*(j-1)+1;
         bigk(rbk:rbk+n1, cbk:cbk+n1) = bigk(rbk:rbk+n1, cbk:cbk+n1) + ke(re:re+n1, ce:ce+n1);
      end
   end
end

% enforce boundary conditions 
% essential bcs assumed homogeneous
for n=1:2*numnod
  if ifix(n)
     for i=1:2*numnod
       force(i) = force(i) - bigk(i,n)*ubar(n);
     end
     bigk(n,:) = zeros(1,numeqns);
     bigk(:,n) = zeros(numeqns,1);
     bigk(n,n) = 1.0;
     force(n) = ubar(n);
  end
end

%solve stiffness equations
disp = force/bigk;

%compute stresses at center of each element
gauss = [0, 0];
ratio = 0;
xc = 0;
yc = 0;
for e=1:numele
   [xcen,ycen,ratioe] = post_process_sigmac(node,x,y,gauss,lambda,mu,sigmac,e,disp);
   if ratioe > ratio
       ratio = ratioe;
       xc = xcen;
       yc = ycen;
   end
end

fprintf("critical point x = [%.3f, %.3f]\n", xc, yc);
end