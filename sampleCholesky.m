function R = sampleCholesky(lambda,mu,CV,Lc,Lx1,Lx2,Ly1,Ly2,Nx,Ny)
%%
%%correlation function
%C = @(tau, Lc) exp(-tau^2/2/Lc^2);
%C = @(tau, Lc) exp(-abs(tau)/Lc);
%C = @(tau, Lc) (1+sqrt(3)*abs(tau)/Lc)*exp(-sqrt(3)*abs(tau)/Lc);
C = @(tau, Lc) (1+sqrt(5)*abs(tau)/Lc+5*tau^2/3/Lc^2)*exp(-sqrt(5)*abs(tau)/Lc);

%%
%%solving shape and scale parameters
disp("solving for shape and scale parameters.")
eta = 0.2;
error = 1;
while error >= 1e-3
    eta = eta + 1e-4;
    error = abs(gamma(1+2/eta) - (CV^2+1) * gamma(1+1/eta) * gamma(1+1/eta));
end
sigma_lambda = lambda / gamma(1+1/eta);
sigma_mu = mu / gamma(1+1/eta);
fprintf("eta = %.4f, sigma_lambda = %.4f, sigma_mu = %.4f\n", eta, sigma_lambda, sigma_mu);

%%
%%mesh
Np = Nx*Ny;
Sx = linspace(Lx1,Lx2,Nx);
Sy = linspace(Ly1,Ly2,Ny);
[X,Y] = meshgrid(Sx,Sy);

%%
%%correlation matrix
disp("setting up correlation matrix...")
pause(.1);
R = eye(Np,Np);
count = 0;
progress = 0;
p_step = 0.1;
fprintf("%2.2f%%\n", progress*100);
for i = 1:Np
    for j = (i+1):Np
        R(i,j) = abs(C(X(i)-X(j),Lc)*C(Y(i)-Y(j),Lc));
        count = count + 1;
    end
    if 2*count/Np/(Np-1) >= (progress + p_step)
        progress = progress + p_step;
        fprintf("%2.2f%%\n", progress*100);
        pause(.1);
    end
end
R = R + R' - diag(diag(R));

%%
%%Cholesky
disp("performing Cholesky factorization...")
pause(.1);
L = chol(R);

%%
%%sampling
disp("sampling...")
pause(.1);
G_lambda = L'*randn(Np,1);
G_mu = L'*randn(Np,1);

%%
%%Weibull transformation
disp("performing Weibull transformation...")
pause(.1);
Z_lambda = wblinv(normcdf(G_lambda,0,1),sigma_lambda,eta);
Z_mu = wblinv(normcdf(G_mu,0,1),sigma_mu,eta);

%%
%%writing lambda
disp("writing lambda to txt...")
pause(.1);
fileID = fopen('lambda.txt', 'w');

fprintf(fileID, 'AXIS X\n');
for i = 1:Nx
    fprintf(fileID, '%.4f ', X(1,i));
end
fprintf(fileID, '\n');

fprintf(fileID, 'AXIS Y\n');
for i = 1:Ny
    fprintf(fileID, '%.4f ', Y(i,1));
end
fprintf(fileID, '\n');

fprintf(fileID, 'DATA\n');
for i = 1:Np
    fprintf(fileID, '%.4f ', Z_lambda(i));
end
fclose(fileID);

%%
%%writing mu
disp("writing mu to txt...")
pause(.1);
fileID = fopen('mu.txt', 'w');

fprintf(fileID, 'AXIS X\n');
for i = 1:Nx
    fprintf(fileID, '%.4f ', X(1,i));
end
fprintf(fileID, '\n');

fprintf(fileID, 'AXIS Y\n');
for i = 1:Ny
    fprintf(fileID, '%.4f ', Y(i,1));
end
fprintf(fileID, '\n');

fprintf(fileID, 'DATA\n');
for i = 1:Np
    fprintf(fileID, '%.4f ', Z_mu(i));
end
fclose(fileID);


%%
%%plot
close all
Z_lambda = reshape(Z_lambda,Ny,Nx);
Z_mu = reshape(Z_mu,Ny,Nx);
figure(1);
surf(X,Y,Z_lambda);
figure(2);
surf(X,Y,Z_mu);

end
