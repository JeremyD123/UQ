N1 = [];
N2 = [];
N3 = [];

lambda_avg = 1.11;
mu_avg = 1.67;
CV = 0.1;
Lc = 25;
L = 100;
H = 100;
Nx = 25;
Ny = 25;
Np = Nx*Ny;
Sx = linspace(0,L,Nx);
Sy = linspace(0,H,Ny);
[Xmesh,Ymesh] = meshgrid(Sx,Sy);
tol = 0.1;

A = 1/CV^2;
B_lambda = lambda_avg*CV^2;
B_mu = mu_avg*CV^2;

[d,v] = KLexpansion(1,Lc,Xmesh,Ymesh,Np,tol);

Nmc = 1000;
nu = length(d);
m = 2*nu;
Y = randn(Nmc,m);

G_lambda = {};
G_mu = {};
lambda = {};
mu = {};

for i = 1:Nmc
    fprintf('i = %d\n', i);
    eta_lambda = Y(i,1:nu)';
    eta_mu = Y(i,nu+1:end)';
    G_lambda{i} = v * (eta_lambda.*sqrt(d));
    G_mu{i} = v * (eta_mu.*sqrt(d));
    
    lambda{i} = gaminv(normcdf(G_lambda{i},0,1),A,B_lambda);
    mu{i} = gaminv(normcdf(G_mu{i},0,1),A,B_mu);
end

G_lambda = cell2mat(G_lambda);
G_mu = cell2mat(G_mu);
lambda = cell2mat(lambda);
mu = cell2mat(mu);

figure
subplot(1,2,1)
[f,x] = ksdensity(G_lambda(:));
plot(x,f,'LineWidth',1.5)
hold on
[f,x] = ksdensity(lambda(:));
plot(x,f,'LineWidth',1.5)
hold off
xlabel('G, \Lambda')
ylabel('P_{G}(x), P_{\Lambda}(x)')
legend('PDF of G_1(x)','PDF of \Lambda(x)')
ax = gca;
ax.FontSize = 24;

subplot(1,2,2)
[f,x] = ksdensity(G_mu(:));
plot(x,f,'LineWidth',1.5)
hold on
[f,x] = ksdensity(mu(:));
plot(x,f,'LineWidth',1.5)
hold off
xlabel('G, M')
ylabel('P_{G}(x), P_{M}(x)')
legend('PDF of G_2(x)','PDF of M(x)')
ax = gca;
ax.FontSize = 24;