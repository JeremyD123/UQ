% MESH GENERATION PROGRAM FOR 2D BEAM IN BENDING

function[x,y,node,numele,numnod] = mesh2d(length,height,Nx,Ny)

% number of elements in each direction
ndivl = Nx;
ndivw = Ny;

numele = ndivw*ndivl;
numnod = (ndivl+1)*(ndivw+1);


% set up nodal coordinates

for i = 1:(ndivl+1)
   for j=1:(ndivw+1)
      x((ndivw+1)*(i-1)+j) = (length/ndivl)*(i-1);
      y((ndivw+1)*(i-1)+j) = height - (height/ndivw)*(j-1);
   end
end

% set up connectivity array

for j=1:ndivl
   for i=1:ndivw
      elemn = (j-1)*ndivw + i;
      nodet(elemn,1) = elemn + (j-1);
      nodet(elemn,2) = nodet(elemn,1) + 1;
      nodet(elemn,3) = nodet(elemn,2) + ndivw + 1;
      nodet(elemn,4) = nodet(elemn,3) -1;
   end
end

node = nodet';