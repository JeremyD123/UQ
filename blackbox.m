function U = blackbox(L,H,Nx,Ny,Y,d,v)
disp('---------------------IN BLACK BOX!---------------------')

%num of KL expansion terms per field
nu = length(d);

%stochastic input for lambda and mu
eta_lambda = Y(1:nu)';
eta_mu = Y(nu+1:end)';

%realization of fields
G_lambda = v * (eta_lambda.*sqrt(d));
G_mu = v * (eta_mu.*sqrt(d));

disp('performing Gamma transformation......')
lambda = gaminv(normcdf(G_lambda,mean(G_lambda),std(G_lambda)),A,B_lambda);
mu = gaminv(normcdf(G_mu,mean(G_lambda),std(G_lambda)),A,B_mu);

U = FEM2D(lambda,mu,L,H,Nx,Ny);

end