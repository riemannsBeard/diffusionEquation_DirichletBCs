clc
clear all
close all

%% Datos

Re = 1;
CFL = 0.5;
N = 128;
L = 1;
dx = L/(N-1);
dt = CFL*Re*dx*dx;
alpha = 0.5*dt/Re;
tf = 5;

x = linspace(0, L, N)';


%% Matrices de diferenciacion

e = ones(N,1);

Dxx = (1/(dx*dx))*spdiags([-e, 2*e, -e], [-1 0 1], N, N);

DxxL = speye(N) + Dxx*alpha;
DxxR = speye(N) - Dxx*alpha;

%% Condiciones de contorno tipo Dirichlet

Uw = 0;
Ue = 0;

bc = zeros(N,1);
bc(1) = Uw;
bc(end) = Ue;
bc = bc/(dx*dx);

% Solucion inicial
x1 = x(x < 0.5);
x2 = x(x >= 0.5);

u = [2*x1; 2*(1 - x2)];

%% Simulación

t = 0;

while (t <= 0.048)
    
    RHS = DxxR*u + bc;
    
    u = DxxL\RHS;
    
    t = t + dt;
    
    u(1) = 0;
    u(end) = 0;
    
%     figure(1)
%     plot(x, u)
%     ylim([0 1])
%     
%     drawnow
%     pause(0.5)
    
end

uTeo = fTeo(x, 0.048, 32);

figure(1)
plot(x, u, x, uTeo)


