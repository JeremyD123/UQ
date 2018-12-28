clear
clc
load('MC_28.mat')
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
% Nmc = 2e4;
% %draw samples, and solve the forward problem Nmc times
% Y = randn(Nmc,m);
% U_mc = [];
% nstars = 0;
% nspaces = 0;
% for i = 1:Nmc
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
%     %U_mc = [U_mc; FEM2D(lambda,mu,L,H,Nx,Ny)];
%     U_mc = [U_mc; calculateStrainEnergy(lambda,mu,Nx*Ny)];
%     disp('solve the forward problem using FEM......')
% end
% 
% %compute convergence
% U_convergence = cumsum(U_mc)'./(1:Nmc);
% U2_convergence = cumsum(U_mc.^2)'./(1:Nmc);
% 
% %plot results
% figure
% subplot(1,3,1);
% plot(U_convergence);
% xlabel('N_{MC}');
% ylabel('E[\Pi]');
% ylim([5.55,5.6]);
% 
% subplot(1,3,2);
% plot(U2_convergence);
% xlabel('N_{MC}');
% ylabel('E[\Pi^2]');
% ylim([30.5,31.5]);
% 
% subplot(1,3,3);
% ksdensity(U_mc);
% xlabel('\pi');
% ylabel('P_\Pi(\pi)')

%%
%%PCE, PSEUDO-PROJECTION
% disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('+                           PCE                           +')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('+                  MONTE CARLO SUMMATION                  +')
% disp('|                                                         |')
% disp('|                                                         |')
% disp('+----+-----+-----+-----+-----+-----+-----+-----+-----+----+')
% %polynomial order
% Q = 3;
% N2 = factorial(m+(0:Q))./factorial(0:Q)/factorial(m);
% N1 = [1,N2(1:end-1)+1];
% %mapping of multi-index
% M_pp = multi_index(m,Q);
% 
% %row = numb of KL terms
% %col = num of RVs, which tells us about the combination of bases
% [row_pp, col_pp] = size(M_pp);
% 
% u_pp = [];
% U_pp = zeros(Nmc,1);
% 
% nstars = 0;
% nspaces = -15;
% progress = 0;
% p_step = 0.02;
% figure
% subplot(1,3,3)
% [f,xi] = ksdensity(U_mc);
% plot(xi,f,'k-','LineWidth',1.5)
% hold on
% for q = 0:Q
%     u_pp_convergence = zeros(1,Nmc);
%     U_pp_convergence = zeros(Nmc,1);
%     for r = N1(q+1):N2(q+1)
%         if r/row_pp >= progress
%             progress = progress + p_step;
%             fprintf(repmat('\b',1,nstars+nspaces+15));
%             nstars = round(r/row_pp*50);
%             nspaces = 50-nstars;
%             fprintf('progress: ||');
%             fprintf(repmat('*',1,nstars));
%             fprintf(repmat('-',1,nspaces));
%             fprintf('||\n');
%         end
%         PSI = 1;
%         for c = 1:col_pp
%             PSI = PSI .* hermiteN(M_pp(r,c),Y(:,c));
%         end
%         %approximate each PCE coefficient
%         u_pp_convergence = u_pp_convergence + (cumsum(U_mc.*PSI)'./(1:Nmc)).^2;
%         u_i = sum(U_mc.*PSI)/Nmc;
%         u_pp = [u_pp, u_i];
%         U_pp = U_pp + u_i*PSI;
%     end
%     subplot(1,3,2)
%     u_pp_convergence = sqrt(u_pp_convergence);
%     plot(u_pp_convergence,'LineWidth',1.5)
%     hold on
%     
%     subplot(1,3,3)
%     [f,xi] = ksdensity(U_pp);
%     plot(xi,f,'LineWidth',1.5)
%     hold on
% end
% subplot(1,3,2)
% hold off
% legend('|\pi^{<0>}|', '|\pi^{<1>}|', '|\pi^{<2>}|', '|\pi^{<3>}|')
% xlabel('N_{MC}')
% ylabel('|\pi^{<\alpha>}|')
% ylim([0,6])
% ax = gca;
% ax.FontSize = 14;
% 
% subplot(1,3,3)
% hold off
% legend('Monte Carlo', 'Q = 0', 'Q = 1', 'Q = 2', 'Q = 3')
% xlabel('\pi')
% ylabel('P_\Pi(\pi)')
% xlim([0,10])
% ax = gca;
% ax.FontSize = 14;
% 
% u_pp_Q = [];
% for i = 0:Q
%     u_pp_Q = [u_pp_Q, norm(u_pp(N1(i+1):N2(i+1)))];
% end
% 
% subplot(1,3,1)
% plot(0:Q,u_pp_Q,'ko-')
% xlabel('|\alpha|')
% ylabel('|\pi^{<\alpha>}|')
% ax = gca;
% ax.FontSize = 14;

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
Q = 3;
N2 = factorial(m+(0:Q))./factorial(0:Q)/factorial(m);
N1 = [1,N2(1:end-1)+1];

%mapping of multi-index
M_sg = multi_index(m,Q);

%row = numb of KL terms
%col = num of RVs, which tells us about the combination of bases
[row_sg, col_sg] = size(M_sg);

%setup sparse grid
[Y_sg,w] = nwspgr('KPN',m,Q+1);
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
    U_colloc = [U_colloc; calculateStrainEnergy(lambda,mu,Nx*Ny)];
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

figure
subplot(1,2,2)
[f,xi] = ksdensity(U_mc);
plot(xi,f,'k-','LineWidth',1.5)
hold on

for q = 0:Q
    U_sg_convergence = zeros(Nmc,1);
    for r = N1(q+1):N2(q+1)
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
    subplot(1,2,2)
    [f,xi] = ksdensity(U_sg);
    plot(xi,f,'LineWidth',1.5)
    hold on
end

%plot results
subplot(1,2,2)
hold off
legend('Monte Carlo', 'Q = 0', 'Q = 1', 'Q = 2', 'Q = 3')
xlabel('\pi')
ylabel('P_\Pi(\pi)')
xlim([0,10])
ax = gca;
ax.FontSize = 14;

u_sg_Q = [];
for i = 0:Q
    u_sg_Q = [u_sg_Q, norm(u_sg(N1(i+1):N2(i+1)))];
end

subplot(1,2,1)
plot(0:Q,u_sg_Q,'ko-')
xlabel('|\alpha|')
ylabel('|\pi^{<\alpha>}|')
ax = gca;
ax.FontSize = 14;