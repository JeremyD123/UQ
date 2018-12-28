clear
clc
load('MC_6.mat')
%%
%%problem setup
% disp('-----+-----+-----+-----+-----+-----+-----+-----+-----+-----')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('+                    SETTING UP PROBLEM                   +')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('+                                                         +')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('-----+-----+-----+-----+-----+-----+-----+-----+-----+-----')
% lambda_avg = 1.11;
% mu_avg = 1.67;
% CV = 0.1;
% Lc = 25;
% L = 100;
% H = 100;
% Nx = 25;
% Ny = 25;
% 
% %solving shape and scale parameters
% disp("solving for shape and scale parameters.")
% A = 1/CV^2;
% B_lambda = lambda_avg*CV^2;
% B_mu = mu_avg*CV^2;
% 
% %mesh
% Np = Nx*Ny;
% Sx = linspace(0,L,Nx);
% Sy = linspace(0,H,Ny);
% [X,Y] = meshgrid(Sx,Sy);
% 
% %KL expansion
% disp("performing KL expansion...")
% tol = 0.1;
% [d,v] = KLexpansion(1,Lc,X,Y,Np,tol);
% 
% %stochastic dimension = num of KL terms per field x num of fields
% nu = length(d);
% m = 2*nu;

%%
%%MONTE CARLO
% disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('+                       MONTE CARLO                       +')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('+                           INIT                          +')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
% 
% %number of samples, Monte Carlo
% Nmc = 1e4;
% %draw samples, and solve the forward problem Nmc times
% Y = randn(Nmc,m);
% U_mc = [];
% nstars = 0;
% nspaces = 0;
% for i = 1:Nmc
%     disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
%     disp('|                                                         |')
%     disp('|                                                         |')
%     disp('+                       MONTE CARLO                       +')
%     disp('|                                                         |')
%     disp('|                                                         |')
%     fprintf('+                      %9d                          +\n',i)
%     disp('|                                                         |')
%     disp('|                                                         |')
%     disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
%     nstars = round(i/Nmc*50);
%     nspaces = 50-nstars;
%     fprintf('progress: ||');
%     fprintf(repmat('*',1,nstars));
%     fprintf(repmat('-',1,nspaces))
%     fprintf('||\n')
%     disp('sampling......') %already sampled :)
%     disp('constructing random fields......')
%     eta_lambda = Y(i,1:nu)';
%     eta_mu = Y(i,nu+1:end)';
%     G_lambda = v * (eta_lambda.*sqrt(d));
%     G_mu = v * (eta_mu.*sqrt(d));
%     
%     disp('performing Gamma transformation......')
%     lambda = gaminv(normcdf(G_lambda,0,1),A,B_lambda);
%     mu = gaminv(normcdf(G_mu,0,1),A,B_mu);
%     fprintf('field statistics:\n')
%     fprintf('lambda: mean = %.4f, std = %.4f\n',mean(lambda),std(lambda));
%     fprintf('    mu: mean = %.4f, std = %.4f\n',mean(mu),std(mu));
%     U_mc = [U_mc; FEM2D(lambda,mu,L,H,Nx,Ny)];
%     disp('solve the forward problem using FEM......')
% end
% 
% %compute convergence
% U_convergence = cumsum(U_mc)'./(1:Nmc);
% U2_convergence = cumsum(U_mc.^2)'./(1:Nmc);
% 
% %plot results
% figure
% subplot(1,2,1);
% plot(U_convergence);
% xlabel('N_{MC}');
% ylabel('E[U]');
% 
% subplot(1,2,2);
% plot(U2_convergence);
% xlabel('N_{MC}');
% ylabel('E[U^2]');

%%
%%PCE, PSEUDO-PROJECTION
disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
disp('|                                                         |')
disp('|                                                         |')
disp('+                           PCE                           +')
disp('|                                                         |')
disp('|                                                         |')
disp('+                  MONTE CARLO SUMMATION                  +')
disp('|                                                         |')
disp('|                                                         |')
disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
%polynomial order
Q = 2;

%mapping of multi-index
M_pp = multi_index(m,Q);

%row = numb of KL terms
%col = num of RVs, which tells us about the combination of bases
[row_pp, col_pp] = size(M_pp);

u_pp = [];
U_pp = zeros(Nmc,1);

nstars = 0;
nspaces = -15;
progress = 0;
p_step = 0.02;
for r = 1:row_pp
    if r/row_pp >= progress
        progress = progress + p_step;
        fprintf(repmat('\b',1,nstars+nspaces+15));
        nstars = round(r/row_pp*50);
        nspaces = 50-nstars;
        fprintf('progress: ||');
        fprintf(repmat('*',1,nstars));
        fprintf(repmat('-',1,nspaces));
        fprintf('||\n');
    end
    PSI = 1;
    for c = 1:col_pp
        PSI = PSI .* hermiteN(M_pp(r,c),Y(:,c));
    end
    %approximate each PCE coefficient
    u_i = sum(U_mc.*PSI)/Nmc;
    u_pp = [u_pp, u_i];
    U_pp = U_pp + u_i*PSI;
end

%%
%%sparse grid
disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
disp('|                                                         |')
disp('|                                                         |')
disp('+                           PCE                           +')
disp('|                        SPARSE GRID                      |')
disp('|                                                         |')
disp('+                                                         +')
disp('|                                                         |')
disp('|                                                         |')
disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
%draw random samples to visualize output
Y = randn(Nmc,m);

%polynomial order
Q = 2;

%mapping of multi-index
M_sg = multi_index(m,Q);

%row = numb of KL terms
%col = num of RVs, which tells us about the combination of bases
[row_sg, col_sg] = size(M_sg);

%setup sparse grid
[Y_sg,w] = nwspgr('KPN',m,2);
[Ncolloc,~] = size(Y_sg);

%evaluate output at each collocation point
U_colloc = [];
for i = 1:Ncolloc
    disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
    disp('|                                                         |')
    disp('|                           PCE                           |')
    disp('+               EVALUATING COLLOCATION POINT              +')
    disp('|                                                         |')
    disp('|                                                         |')
    fprintf('+                        %5d                            +\n',i)
    fprintf('|                    out of %5d                         |\n',Ncolloc)
    disp('|                                                         |')
    disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
    nstars = round(i/Ncolloc*50);
    nspaces = 50-nstars;
    fprintf('progress: ||');
    fprintf(repmat('*',1,nstars));
    fprintf(repmat('-',1,nspaces))
    fprintf('||\n')
    disp('sampling......') %already sampled :)
    disp('constructing random fields......')
    eta_lambda = Y_sg(i,1:nu)';
    eta_mu = Y_sg(i,nu+1:end)';
    G_lambda = v * (eta_lambda.*sqrt(d));
    G_mu = v * (eta_mu.*sqrt(d));
    
    disp('performing Gamma transformation......')
    lambda = gaminv(normcdf(G_lambda,0,1),A,B_lambda);
    mu = gaminv(normcdf(G_mu,0,1),A,B_mu);
    fprintf('field statistics:\n')
    fprintf('lambda: mean = %.4f, std = %.4f\n',mean(lambda),std(lambda));
    fprintf('    mu: mean = %.4f, std = %.4f\n',mean(mu),std(mu));
    U_colloc = [U_colloc; FEM2D(lambda,mu,L,H,Nx,Ny)];
    disp('solve the forward problem using FEM......')
end

%evaluate PCE coefficients using collocation
disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
disp('|                                                         |')
disp('|                                                         |')
disp('+                           PCE                           +')
disp('|                        SPARSE GRID                      |')
disp('|                                                         |')
disp('+                  EVALUATING COEFFICIENTS                +')
disp('|                                                         |')
disp('|                                                         |')
disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
u_sg = [];
U_sg = zeros(Nmc,1);

nstars = 0;
nspaces = -15;
progress = 0;
p_step = 0.02;
for r = 1:row_sg
    if r/row_sg >= progress
        progress = progress + p_step;
        fprintf(repmat('\b',1,nstars+nspaces+15));
        nstars = round(r/row_sg*50);
        nspaces = 50-nstars;
        fprintf('progress: ||');
        fprintf(repmat('*',1,nstars));
        fprintf(repmat('-',1,nspaces));
        fprintf('||\n');
    end
    PSI = 1;
    PSI_mc = 1;
    for c = 1:col_sg
        PSI = PSI .* hermiteN(M_sg(r,c),Y_sg(:,c));
        PSI_mc = PSI_mc .* hermiteN(M_sg(r,c),Y(:,c));
    end
    %approximate each PCE coefficient
    u_i = w'*(U_colloc.*PSI);
    u_sg = [u_sg, u_i];
    U_sg = U_sg + u_i*PSI_mc;
end

%plot results
figure
subplot(1,3,1);
plot(0:row_pp-1,u_pp,'LineWidth',3);
xlabel('i')
ylabel('u_i')
title('convergence of coeff evaluated using MC summation')

subplot(1,3,2);
plot(0:row_sg-1,u_sg,'LineWidth',3);
xlabel('i')
ylabel('u_i')
title('convergence of coeff using sparse grid quadrature')

subplot(1,3,3);
[P_U_mc,U_mc] = ksdensity(U_mc);
[P_U_pp,U_pp] = ksdensity(U_pp);
[P_U_sg,U_sg] = ksdensity(U_sg);
plot(U_mc,P_U_mc,'k','LineWidth',2)
hold on
plot(U_pp,P_U_pp,'b--','LineWidth',2)
plot(U_sg,P_U_sg,'r--','LineWidth',2)
hold off
xlabel('U')
ylabel('P_U(Y)')
legend('Monte Carlo','PCE-Monte Carlo','PCE-sparse grid quadrature')