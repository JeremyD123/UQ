clear
clc

%problem parameters
b = 0.5;
h = 0.3;
L = 5;
P = 0.01;

%shape parameter
A = [25,100];
%scale parameter
B = [1200,0.05];

m = length(A);
Q = 3;
M = multi_index(m,Q);
[row, col] = size(M);

%%
%%MONTE CARLO

%number of samples, Monte Carlo
Nmc0 = 1e6;
%draw random samples Nmc times
Y = [];
for c = 1:col
    Y = [Y,gamrnd(A(c),B(c),Nmc0,1)];
end
%solve the problem for each sample
U = @(Y) P.*Y(:,2).^3/4/b/h^3./Y(:,1);
U_mc = U(Y);
U_convergence = cumsum(U_mc)'./(1:Nmc0);
figure(1);
subplot(1,2,1);
plot(U_convergence);
subplot(1,2,2);
ksdensity(U_mc);
hold on

%%
%%PCE, PSEUDO-PROJECTION

% number of samples, Monte Carlo
Nmc = 1e6;
Npce = factorial(m+Q)/factorial(m)/factorial(Q);
u_pp = [];

%draw random samples Nmc times
X = [];
for c = 1:col
    X = [X,gamrnd(A(c),1,Nmc,1)];
end
U_pp = zeros(Nmc,1);

for r = 1:row
    PSI = 1;
    for c = 1:col
        PSI = PSI .* laguerreN(M(r,c),A(c)-1,X(:,c));
    end
    fprintf("finding u_{%d}.............\n",r);
    %approximate each PCE coefficient
    u_i = sum(U(X*diag(B)).*PSI)/Nmc;
    u_pp = [u_pp, u_i];
    U_pp = U_pp + u_i*PSI;
end

%plot convergence of PCE coefficients
ksdensity(U_pp);
legend('Monte Carlo','pseudo projection');
hold off

%%
%%PCE, SPARSE GRID
disp('====================PCE, SPARSE GRID======================')

