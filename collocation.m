%problem parameters
b = 0.5;
h = 0.3;
L = 5;
P = 0.01;

%shape parameter
A = 25;
%scale parameter
B = 1200;

U = @(X) P*L^3/4/b/h^3/B./X;
psi = @(i,X) sqrt(factorial(i)*gamma(A)/gamma(i+A))*laguerreL(i,A-1,X)/gamma(A);

u = [];
for i = 0:8
    [z,w] = gen_laguerre_rule(20,A-1,0,1,'qr');
    u = [u, sum(U(z).*w.*psi(i,z))];
end