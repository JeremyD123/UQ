%%
%%term minimization

N1 = [];
N2 = [];
N3 = [];

L = 50;
H = 50;
Nx = 25;
Ny = 25;
Np = Nx*Ny;
Sx = linspace(0,L,Nx);
Sy = linspace(0,H,Ny);
[Xmesh,Ymesh] = meshgrid(Sx,Sy);

tol = 0.05;

%square exponential correlation function
disp("-------------------------------------------")
disp("square exponential")
for Lc = 1:50
    fprintf("Lc = %.1f\n", Lc);
    [d,v] = KLexpansion(1,Lc,Xmesh,Ymesh,Np,tol);
    N1 = [N1, length(d)];
end

%matern class v = 3/2
disp("-------------------------------------------")
disp("Matern class \nu = 3/2")
for Lc = 1:50
    fprintf("Lc = %.1f\n", Lc);
    [d,v] = KLexpansion(2,Lc,Xmesh,Ymesh,Np,tol);
    N2 = [N2, length(d)];
end

%matern class v = 5/2
disp("-------------------------------------------")
disp("Matern class \nu = 5/2")
for Lc = 1:50
    fprintf("Lc = %.1f\n", Lc);
    [d,v] = KLexpansion(3,Lc,Xmesh,Ymesh,Np,tol);
    N3 = [N3, length(d)];
end

%plot convergence
figure
subplot(1,2,1)
plot(1:25,N1(1:25),'k.-','LineWidth',1)
hold on
plot(1:25,N2(1:25),'r.-','LineWidth',1)
plot(1:25,N3(1:25),'b.-','LineWidth',1)
hold off
xlabel('Lc')
ylabel('q')
%title('number of expansion terms for 5% error')
legend('Square Exponential','Matern class \nu = 3/2','Matern class \nu = 5/2')
ax = gca;
ax.FontSize = 14;

%%
%%error convergence
e1 = [];
e2 = [];
e3 = [];

L = 50;
H = 50;
Nx = 25;
Ny = 25;
Np = Nx*Ny;
Sx = linspace(0,L,Nx);
Sy = linspace(0,H,Ny);
[Xmesh,Ymesh] = meshgrid(Sx,Sy);

%square exponential correlation function
disp("-------------------------------------------")
disp("square exponential")
for q = 1:100
    fprintf("q = %.1f\n", q);
    [d,v,error] = KLexpansion2(1,5,Xmesh,Ymesh,Np,q);
    e1 = [e1, error];
end

%matern class v = 3/2
disp("-------------------------------------------")
disp("Matern class \nu = 3/2")
for q = 1:100
    fprintf("q = %.1f\n", q);
    [d,v,error] = KLexpansion2(2,5,Xmesh,Ymesh,Np,q);
    e2 = [e2, error];
end

%matern class v = 5/2
disp("-------------------------------------------")
disp("Matern class \nu = 5/2")
for q = 1:100
    fprintf("q = %.1f\n", q);
    [d,v,error] = KLexpansion2(3,5,Xmesh,Ymesh,Np,q);
    e3 = [e3, error];
end

%plot convergence
subplot(1,2,2)
plot(1:100,e1,'k.-','LineWidth',1)
hold on
plot(1:100,e2,'r.-','LineWidth',1)
plot(1:100,e3,'b.-','LineWidth',1)
hold off
xlabel('q')
ylabel('\epsilon(q)')
%title('number of expansion terms for 5% error')
legend('Square Exponential','Matern class \nu = 3/2','Matern class \nu = 5/2')
ax = gca;
ax.FontSize = 14;