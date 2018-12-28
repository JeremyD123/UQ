clear
clc

lambda_avg = 1.11;
mu_avg = 1.67;
CV = 0.1;
Lc = 25;
L = 100;
H = 100;
Nx = 25;
Ny = 25;

%solving shape and scale parameters
A = 1/CV^2;
B_lambda = lambda_avg*CV^2;
B_mu = mu_avg*CV^2;

%mesh
Np = Nx*Ny;
Sx = linspace(0,L,Nx);
Sy = linspace(0,H,Ny);
[Xmesh,Ymesh] = meshgrid(Sx,Sy);

%KL expansion
tol = 0.1;

[d,v] = KLexpansion(1,Lc,Xmesh,Ymesh,Np,tol);
nu = length(d);
m = 2 * nu;

Y = randn(2,m);

eta_lambda = Y(1,1:nu)';
eta_mu = Y(2,nu+1:end)';
G_lambda = v * (eta_lambda.*sqrt(d));
G_mu = v * (eta_mu.*sqrt(d));

lambda = gaminv(normcdf(G_lambda,0,1),A,B_lambda);
mu = gaminv(normcdf(G_mu,0,1),A,B_mu);

G_lambda = reshape(G_lambda,Ny,Nx);
lambda = reshape(lambda,Ny,Nx);

G_mu = reshape(G_mu,Ny,Nx);
mu = reshape(mu,Ny,Nx);

figure
surf(Xmesh,Ymesh,G_lambda)
ax = gca;
ax.FontSize = 72;
figure
surf(Xmesh,Ymesh,lambda)
ax = gca;
ax.FontSize = 72;
figure
surf(Xmesh,Ymesh,G_mu)
ax = gca;
ax.FontSize = 72;
figure
surf(Xmesh,Ymesh,mu)
ax = gca;
ax.FontSize = 72;