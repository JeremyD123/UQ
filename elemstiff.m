function [ke] = elemstiff(node,x,y,gauss,lambda,mu,e)

% 2d QUAD element stiffness routine

ke = zeros(8,8);
one = ones(1,4);
psiJ = [-1, +1, +1, -1]; etaJ = [-1, -1, +1, +1];

% plane stress D matrix
young = mu(e)*(3*lambda(e)+2*mu(e))/(lambda(e)+mu(e));
pr = lambda(e)/2/(lambda(e)+mu(e));
fac = young/(1 - pr^2);
D = fac*[1.0, pr, 0;
         pr, 1.0, 0.0;
         0, 0, (1.-pr)/2 ];
      
% get coordinates of element nodes 
for j=1:4
   je = node(j,e); xe(j) = x(je); ye(j) = y(je);
end

% compute element stiffness
% loop over gauss points in eta
for i=1:2
   % loop over gauss points in psi
   for j=1:2
      eta = gauss(i);  psi = gauss(j);
      % compute derivatives of shape functions in reference coordinates
      NJpsi = 0.25*psiJ.*(one + eta*etaJ);
      NJeta = 0.25*etaJ.*(one + psi*psiJ);
      % compute derivatives of x and y wrt psi and eta
      xpsi = NJpsi*xe'; ypsi = NJpsi*ye'; xeta = NJeta*xe';  yeta = NJeta*ye';
      Jinv = [yeta, -xeta; -ypsi, xpsi]';
      jcob = xpsi*yeta - xeta*ypsi;
      % compute derivatives of shape functions in element coordinates
      NJdpsieta = [NJpsi; NJeta];
      NJdxy = Jinv*NJdpsieta;
      % assemble B matrix
      BJ = zeros(3,8);
      BJ(1,1:2:7) = NJdxy(1,1:4);  BJ(2,2:2:8) = NJdxy(2,1:4);
      BJ(3,1:2:7) = NJdxy(2,1:4);  BJ(3,2:2:8) = NJdxy(1,1:4);
      % assemble ke
      ke = ke + BJ'*D*BJ/jcob;
   end
end
