function Eelastic = FEM2D(lambda,mu,length,height,Nx,Ny)
%INPUT DATA
[x,y,node,numele,numnod] = mesh2d(length,height,Nx,Ny);

%MATERIAL CONSTANTS
%already passed in
%a field, possibly random, of lambda
%a field, possibly random, of mu

%FORCE AND DISPLACEMENT BC'S
g = 1.0;
[force,ifix,ubar] = applybcs(x,y,numnod,length,height,g);

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
stress = zeros(numele,6);
strain = zeros(numele,6);
Eelastic = 0;
for e=1:numele
   [stresse,straine,Eelastice] = post_process(node,x,y,gauss,lambda,mu,e,disp);
   stress(e,1:6) = stresse;
   strain(e,1:6) = straine;
   Eelastic = Eelastic + Eelastice;
end

% %visualize stress and strain
% figure(1)
% subplot(1,2,1)
% scatter(stress(:,2),stress(:,3),20,(stress(:,4)+stress(:,5))/2,'filled')
% colorbar
% subplot(1,2,2)
% scatter(strain(:,2),strain(:,3),20,(strain(:,4)+strain(:,5))/2,'filled')
% colorbar

fprintf("total elastic energy = %.4f\n", Eelastic);
end