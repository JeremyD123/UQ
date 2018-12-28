%problem parameters
b = 0.5;
h = 0.3;
L = 5;
P = 0.01;

U = @(Y) P*L.^3/4/b/h^3./Y;
%shape parameter
A = 25;
%scale parameter
B = 1200;

%%
%%MONTE CARLO
disp('======================MONTE CARLO========================')

%number of samples, Monte Carlo
Nmc0 = 1e6;
%draw random samples Nmc times
disp('sampling..........')
Y = gamrnd(A,B,Nmc0,1);
%solve the problem for each sample
U_mc = P*L^3./(4*b*h^3*Y);
%compute convergence
U_convergence = cumsum(U_mc)'./(1:Nmc0);
U2_convergence = cumsum(U_mc.^2)'./(1:Nmc0);

%plot results
figure(1);
subplot(1,2,1);
plot(U_convergence);
xlabel('N_{MC}');
ylabel('E[U]');

subplot(1,2,2);
plot(U2_convergence);
xlabel('N_{MC}');
ylabel('E[U^2]');

%%
%%PCE, PSEUDO-PROJECTION
disp('================PCE, PSEUDO PROJECTION===================')

% number of samples, Monte Carlo
Nmc = 1e6;
Npce = 9;
u_sp = [];

%draw random samples Nmc times
disp('sampling..........')
X = gamrnd(A,1,Nmc,1);
U_sp = zeros(Nmc,1);

for i = 0:Npce-1
    fprintf("finding u%d.............\n",i);
    %approximate each PCE coefficient
    u_i = sum(U(B*X).*laguerreN(i,A-1,X))/Nmc;
    u_sp = [u_sp, u_i];
    U_sp = U_sp + u_i*laguerreN(i,A-1,X);
end

%plot results
figure(2);
subplot(1,2,1);
plot(0:Npce-1,u_sp,'o-');
hold on

subplot(1,2,2);
ksdensity(U_mc);
hold on
ksdensity(U_sp);

%%
%%collocation
disp('====================PCE, QUADRATURE======================')

% number of samples, Monte Carlo
Nmc = 1e6;
Npce = 9;
u_colloc = [];

%draw random samples Nmc times
X = gamrnd(A,1,Nmc,1);
U_colloc = zeros(Nmc,1);

for i = 0:Npce-1
    fprintf("finding u%d.............\n",i);
    [z,w] = gen_laguerre_rule(20,A-1,0,1,'qr');
    u_i = sum(U(B*z).*w.*laguerreN(i,A-1,z))/gamma(A);
    u_colloc = [u_colloc, u_i];
    U_colloc = U_colloc + u_i*laguerreN(i,A-1,X);
end

%plot results
subplot(1,2,1);
plot(0:Npce-1,u_colloc,'o-');

subplot(1,2,2);
ksdensity(U_colloc);

%%
%%sparse grid
disp('====================PCE, SPARSE GRID======================')

% number of samples, Monte Carlo
Nmc = 1e6;
Npce = 9;
u_sg = [];

%draw random samples Nmc times
X = gamrnd(A,1,Nmc,1);
U_sg = zeros(Nmc,1);
options = spset('MinDepth',6,'MaxDepth',6,'GridType','Gauss-Patterson','SparseIndices','on');

for i = 0:Npce-1
    fprintf("finding u%d.............\n",i);
    f = @(X) (U(B*X).*laguerreN(i,A-1,X).*X.^(A-1).*exp(-X) + U(B./X).*laguerreN(i,A-1,1./X).*X.^(-A-1).*exp(-1./X))/gamma(A);
    z = spvals(f,1,[],options);
    u_i = spquad(z);
    u_sg = [u_sg, u_i];
    U_sg = U_sg + u_i*laguerreN(i,A-1,X);
end

%plot results
subplot(1,2,1);
plot(0:Npce-1,u_sg,'o-');
ylabel('u_i');
xlabel('i');
legend('pseudo-projection','collocation','sparse grid');

subplot(1,2,2);
ksdensity(U_sg);
ylabel('P_Y');
xlabel('U(Y)');
xlim([0,2e-3]);
legend('Monte Carlo','pseudo-projection','collocation','sparse grid');