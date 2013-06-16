% Heat equation on domain (0,1) using semidiscrete method
%
% Dale Roberts <dale.o.roberts@gmail.com>

clear all;

%% Parameters

nx=10;

%% Create mesh

dx=1/(nx-1);
x=0:dx:1;

%% Initial condition

u=sin(pi*x)';

%% 1D Laplacian

e=ones(nx,1);
A=1/(dx^2)*spdiags([e -2*e e],-1:1,nx,nx);

%% Setup problem and solve

global A

options=odeset('AbsTol',1e-10,'RelTol',1e-10);
[t U]=ode15s('Au',[0 1],u,options);

%% Plot solution

% style
set(gcf, 'Renderer', 'opengl');
shading faceted
colormap Jet

surf(x,t,U,'EdgeColor','None');
view(140, 30);
xlabel('space x');
ylabel('time t');
zlabel('u');
