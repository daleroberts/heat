% Solve heat equation using Crank-Nicholson
%
% Dale Roberts <dale.o.roberts@gmail.com>

clear all;

%% Parameters

nx=30;
nt=40;

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

r=dt/(dx^2);
T=speye(nx)-r*A;

U=zeros(nx+2,nt);
U(2:end-1,1) = u;

for j=2:nt
    b = u+r*A*u;
    u = T\b;
    U(2:end-1,j) = u;
end

%% Plot solution

% style
set(gcf, 'Renderer', 'opengl');
shading faceted
colormap Jet

surf(x,t,U(2:end-1,:)','EdgeColor','None');
view(140, 30);
xlabel('space x');
ylabel('time t');
zlabel('u');
