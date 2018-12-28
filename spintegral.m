%problem parameters
b = 0.5;
h = 0.3;
L = 5;
P = 0.01;
A = 25;
B = 1200;

d = 1;
Npce = 9;
U = @(Y) P*L^3/4/b/h^3./Y;
u_pce = [];

for i = 0:Npce-1
    f = @(X) (U(B*X).*laguerreN(i,A-1,X).*X.^(A-1).*exp(-X) + U(B./X).*laguerreN(i,A-1,1./X).*X.^(-A-1).*exp(-1./X))/gamma(A);
    [x,w] = spquad(1,3,[0;1]);
    q = w'*f(x);
    u_pce = [u_pce, q];
end

plot(0:Npce-1,u_pce,'o-')