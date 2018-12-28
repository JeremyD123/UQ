clear
clc

lambda_avg = 1.11;
mu_avg = 1.67;
CV = 0.1;
Lc = 5;
L = 100;
H = 100;
Nx = 50;
Ny = 50;

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
tol = 0.01;

t = {'Exponential', 'Square Exponential', 'Matern \nu = 3/2', 'Matern \nu = 5/2'};
figure;

for i = 0:3
    [d,v] = KLexpansion(i,Lc,Xmesh,Ymesh,Np,tol);
    
    %stochastic dimension = num of KL terms per field x num of fields
    nu = length(d);
    
    Y = randn(1,nu);
    eta = Y(1,:)';
    G = v * (eta.*sqrt(d));
    
    G = reshape(G,Ny,Nx);
    
    subplot(2,2,i+1)
    surf(Xmesh,Ymesh,G);
    title(t{i+1})
end