% Explicit finite difference for the heat equation on domain (0,1)
%
% Dale Roberts <dale.o.roberts@gmail.com>

clear all;

%% Parameters

nx=10;
nt=200;

%% Create mesh

dx=1/(nx-1);
dt=1/(nt-1);

x=0:dx:1;
t=0:dt:1;

%% Initial condition

u=sin(pi*x)';

%% 1D Laplacian

e=ones(nx,1);
A=spdiags([e -2*e e],-1:1,nx,nx);

%% Setup problem and solve

U=zeros(nx,nt);
U(:,1)=u;

mu=dt/(dx^2) % mu < 0.5?
A=mu*A;
for j=2:nt
    U(:,j)=U(:,j-1)+A*U(:,j-1);
end

%% Plot solution

% style
set(gcf, 'Renderer', 'opengl');
shading faceted
colormap Jet

surf(x,t,U','EdgeColor','None');
view(140, 30);
xlabel('space x');
ylabel('time t');
zlabel('u');
