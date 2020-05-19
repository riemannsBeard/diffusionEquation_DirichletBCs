clc
clear all
close all

%% Datos

Re = 1;
Nx = 256;
Ny = 128;
Lx = 2;
Ly = 1;
tf = 0.05;
CFL = 5;
Qc = 10;

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);

X = linspace(0, Lx, Nx);
Y = linspace(0, Ly, Ny);

dt = min(Re*dx*dx/CFL, Re*dy*dy/CFL);

[x, y] = meshgrid(X, Y);

%% Matrices de diferenciacion

ex = ones(Nx,1);
Dxx = (1/dx^2)*spdiags([ex -2*ex ex], [-1 0 1], Nx, Nx);

ey = ones(Ny,1);
Dyy = (1/dy^2)*spdiags([ey -2*ey ey], [-1 0 1], Ny, Ny);
 
L = kron(Dxx, speye(Ny)) + kron(speye(Nx), Dyy);

u = 0*speye(Ny,Nx);

A = speye(Ny*Nx) - 0.5*dt*L/Re;
B = speye(Ny*Nx) + 0.5*dt*L/Re;

%% Condiciones de contorno tipo Dirichlet

bc = 0*speye(Ny,Nx);

bc(:,1) = Qc; % E
bc(:,Nx) = Qc; % W
bc(1,:) = Qc; % N
bc(Ny,:) = Qc; % S

% Solucion inicial con condiciones de contorno
u = u + bc;

%% Simulacion

% Pinta estado inicial
figure (1),
pcolor(x, y, u)
colormap jet
shading interp;
caxis([0 10])
daspect([1/Lx 1/Ly 1])
colorbar
drawnow

% Redistribuye las matrices
u = reshape(u, Nx*Ny, 1);
bc = reshape(bc, Nx*Ny, 1);

% writerObj = VideoWriter('out.avi'); % Name it.
% writerObj.FrameRate = 1; % How many frames per second.
% open(writerObj);

t = 0; j = 1;

% aviobj = VideoWriter('out.avi', 'Uncompressed avi');
% aviobj.FrameRate = 5;
% open(aviobj);

while t<tf
       
    RHS = B*u + bc/CFL;
    
    u = A\RHS;
    
    t = t + dt;
    
    if rem(j, 10) == 0
        
        u = reshape(u, Ny, Nx);
        
        figure (1),
        %contourf(x, y, u, 8)
        pcolor(x, y, u)
        colormap jet
        shading interp;
        caxis([0 10])
        colorbar
        daspect([1 1 1])
        drawnow
        
%         F = getframe(gcf);
%         writeVideo(aviobj, F);

        u = reshape(u, Ny*Nx, 1);

    end
    
    j = j+1;
end

% close(aviobj);

%% Video


