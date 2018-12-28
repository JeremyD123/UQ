function [Z_lambda,Z_mu,nu] = sampleKL(lambda,mu,CV,Lc,Lx1,Lx2,Ly1,Ly2,Nx,Ny,varargin)
%%
%%solving shape and scale parameters
disp("solving for shape and scale parameters.")
A = 1/CV^2;
B_lambda = lambda*CV^2;
B_mu = mu*CV^2;

%%
%%mesh
Np = Nx*Ny;
Sx = linspace(Lx1,Lx2,Nx);
Sy = linspace(Ly1,Ly2,Ny);
[X,Y] = meshgrid(Sx,Sy);

%%
%%KL expansion
disp("performing KL expansion...")
[d,v] = KLexpansion(1,Lc,X,Y,Np);

%%
%%sampling
if length(varargin) ~= 2
    disp("sampling......")
    nu = length(d);
    eta_lambda = randn(nu,1);
    eta_mu = randn(nu,1);
end
G_lambda = v * (eta_lambda.*sqrt(d));
G_mu = v * (eta_mu.*sqrt(d));
nu = length(d);

%%
%%inverse gamma transformation
disp("performing Gamma transformation...")
Z_lambda = gaminv(normcdf(G_lambda,0,1),A,B_lambda);
Z_mu = gaminv(normcdf(G_mu,0,1),A,B_lambda);

%%
%%plot results
% Z_lambda_reshape = reshape(Z_lambda,Ny,Nx);
% Z_mu_reshape = reshape(Z_mu,Ny,Nx);
% figure(1)
% subplot(1,2,1)
% surf(X,Y,Z_lambda_reshape)
% subplot(1,2,2)
% surf(X,Y,Z_mu_reshape)


end
