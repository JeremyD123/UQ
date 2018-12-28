Lc = 5;
C = @(x, y) exp(-(x.^2+y.^2)./2./Lc^2);
lambda = 500;
behavior = 1;
a = -50;
b = 50;
RHS = @(x) x.*0;
sol = Fie(lambda,a,b,behavior,C,RHS,1e-6,1e-6);
