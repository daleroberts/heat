% ADI finite-difference scheme for the 2D heat equation
%
% Dale Roberts <dale.o.roberts@gmail.com>

close all;
clear all;

%% Set mesh size

xmin=0; xmax=1;
ymin=0; ymax=1;

%% Create uniform mesh of the domain

Nx=9; Ny=39;

d.x=linspace(xmin,xmax,Nx);
d.y=linspace(ymin,ymax,Ny);

%% Setup the domain

[X,Y]=meshgrid(d.x,d.y);

% Identify interior
d.interior=(X>0);
d.interior(1,:)=false;
d.interior(end,:)=false;
d.interior(:,1)=false;
d.interior(:,end)=false;

% Index unknowns and setup lookup table
r=0;
d.rr=zeros(Nx,Ny);
for j=1:Ny
    for i=1:Nx
        if d.interior(i,j)
            r=r+1;
            d.ii(r)=i;
            d.jj(r)=j;
            d.rr(i,j)=r;
        end
    end
end

% Save number of unknowns
d.NU=r;

%% Display the mesh lines

figure('Name','Mesh');
hold all;
axis([xmin xmax ymin ymax]);

% draw vertical lines
ymin=d.y(1); ymax=d.y(end);
for i=1:length(d.x)
    line([d.x(i) d.x(i)],[ymin ymax],'color','k');
end

% draw horizontal lines
xmin=d.x(1); xmax=d.x(end);
for i=1:length(d.y)
    line([xmin xmax],[d.y(i) d.y(i)],'color','k');
end

% draw labels
xlabel('x');
ylabel('y');

%% Display the mesh points

% Mesh points
[X,Y]=meshgrid(d.x,d.y);

% draw interior mesh points
plot(X(d.interior), Y(d.interior),'.','Color','black','MarkerSize',6.0);

% draw edge mesh points
plot(X(~d.interior),Y(~d.interior),'.','Color','blue','MarkerSize',6.0);

%% Annotate the mesh points

for r=1:d.NU
    % get mesh indices
    i=d.ii(r);
    j=d.jj(r);
    text(d.x(i),d.y(j),sprintf('%i',r),'FontSize',9,'HorizontalAlignment','center', ...
        'EdgeColor','black','BackgroundColor','white','VerticalAlignment','middle');
end

%% Finish mesh display

hold off;

%% Generate the x-direction operator

Ax=spalloc(d.NU,d.NU,3*d.NU);

% Second-order central difference in x coordinate
for r=1:d.NU
    % get mesh indices
    i=d.ii(r);
    j=d.jj(r);
    % calculate mesh spacing
    h1=(d.x(i)-d.x(i-1));
    h2=(d.x(i+1)-d.x(i));
    % central element
    Ax(r,r)=-2/(h1*h2);
    % west element
    if d.interior(i-1,j)
        Ax(r,d.rr(i-1,j))=2/(h1*(h1+h2));
    end
    % east element
    if d.interior(i+1,j)
        Ax(r,d.rr(i+1,j))=2/(h2*(h1+h2));
    end
end

% figure('Name','Ax');
% spy(Ax)

%% Generate the y-direction operator

Ay=spalloc(d.NU,d.NU,3*d.NU);

% Second-order central difference in y coordinate
for r=1:d.NU
    % get mesh indices
    i=d.ii(r);
    j=d.jj(r);
    % calculate mesh spacing
    h1=(d.y(j)-d.y(j-1));
    h2=(d.y(j+1)-d.y(j));
    % central element
    Ay(r,r)=-2/(h1*h2);
    % north element
    if d.interior(i,j+1)
        Ay(r,d.rr(i,j+1))=2/(h1*(h1+h2));
    end
    % south element
    if d.interior(i,j-1)
        Ay(r,d.rr(i,j-1))=2/(h2*(h1+h2));
    end
end

% figure('Name','Ay');
% spy(Ay);

%% Set A matrix

A=Ax+Ay;

% figure('Name','A');
% spy(A);

%% Set initial condition

U=zeros(d.NU);

for r=1:d.NU
    % get mesh indices
    i=d.ii(r);
    j=d.jj(r);
    % set initial condition
    U(r)=100.0*sin(pi*d.x(i))*sin(pi*d.y(j));
end

UU=zeros(length(d.x),length(d.y));
for r=1:d.NU
    % get mesh indices
    i=d.ii(r);
    j=d.jj(r);
    % set value
    UU(i,j)=U(r);
end

% figure('Name','Initial condition');
% surf(d.x,d.y,UU)

%% Perform Crank-Nicolson scheme

Nt=10;
dt=0.01;

T=(eye(d.NU)-0.5*dt*A);

for k=1:Nt
    b=U+0.5*dt*A*U;
    U=T\b;
end

t=k*dt;

%% Plot solution

UU=zeros(length(d.x),length(d.y));
for r=1:d.NU
    % get mesh indices
    i=d.ii(r);
    j=d.jj(r);
    % set value
    UU(i,j)=U(r);
end


fprintf('t=%f\n',t);

figure('Name','Solution');
surf(d.x,d.y,UU)

%% Compare against exact solution

u = @(x,y,t) 100*exp(-2*pi^2*t)*sin(pi*x)*sin(pi*y);

VV=zeros(length(d.x),length(d.y));
for r=1:d.NU
    % get mesh indices
    i=d.ii(r);
    j=d.jj(r);
    % set value
    VV(i,j)=abs(U(r)-u(d.x(i),d.y(j),t));
end

figure('Name','Error');
surf(d.x,d.y,VV)

% draw labels
xlabel('x');
ylabel('y');
