clear all
close all

%space variables (ENTER)
Nx=50;           %number of columns
Ny=50;           %number of rows 

%discretization variables (ENTER)
dx = 1;          % x-grid size
dy = 1;          % y-grid size


%fluid density & viscosity (ENTER)
density = 1000;
viscosity = 1;
ratio = density/viscosity;       %ratio of density and dynamic viscosity 

%size and dimension of pressure and velocity (AUTO)
p = zeros(Ny,Nx);      %pressure
u = zeros(Ny+1,Nx+1)+0.1;  %x-velocity  
v = zeros(Ny+1,Nx+2);  %y-velocity
residual = zeros(Ny*Nx,1); %residuals from continuity
dp = zeros(Ny*Nx,1);  %changes in pressures

%initial conditions (ENTER)
u = zeros(Ny+1,Nx+1)+0.1;   %constant value of 0.1 (initializes velocity field)

%temporary variables (AUTO)
u1=u;  
v1=v;
dp1=zeros(Nx,Ny);
residual1=zeros(Nx,Ny);

%timestep value, relaxation factor, number of iterations (ENTER)
dt=1;             
relaxation_factor=0.5;            
total_iterations=1000;            
residual_max = zeros(total_iterations,1);

%check CFL criteria (CHECK!)
CFL_x = max(max(u))*dt/dx;
CFL_y = max(max(u))*dt/dy;

%calculate sparse matrix (AUTO)
J_a = 2*(1/dx^2+1/dy^2);
J_b = -1/dy^2;
J_c = -1/dx^2;

J=spalloc(Nx*Ny,Nx*Ny,(Nx-2)*(Ny-2)*4+Nx*Ny);

for i=1:Nx*Ny-1
    J(i,i+1)=J_b/J_a;
end

for i=2:Nx*Ny
    J(i,i-1)=J_b/J_a;
end

for i=1:Nx*Ny-Nx
    J(i,i+Nx)=J_c/J_a;
end

for i=1:Nx*Ny-Nx
    J(i+Nx,i)=J_c/J_a;
end

for i=1:Ny-1
J(i*Nx+1,i*Nx)=0;
J(i*Nx,i*Nx+1)=0;
end

for i=1:Nx*Ny
    J(i,i)=-sum(J(i,:));
end

%spy(J)   %for checking sparse matrix

    
%calculate velocity field using Navier-Stokes equations    
for j=2:Ny
    for i=2:Nx
        
        a1=-((u(j,i+1))^2-(u(j,i-1))^2)/(2*dx)  -  (u(j+1,i)*(v(j+1,i+1)+v(j+1,i))-u(j-1,i)*(v(j,i+1)+v(j,i)))/(4*dy);
        a3=(u(j,i+1)-2*u(j,i)+u(j,i-1))/(dx^2);            %solving N-S Eq. for u velocity 
        a4=(u(j+1,i)-2*u(j,i)+u(j-1,i))/(dy^2);
          
        A=a1+(a3+a4)/ratio;
      
        u1(j,i)=u(j,i)+dt*(A-(p(j,i)-p(j,i-1))/dx);
        
    end
end

for j=2:Ny
    for i=2:Nx+1
        
        b1=-((((v(j+1,i))^2-(v(j-1,i))^2)/(2*dy))  -  ((v(j,i+1)*(u(j-1,i)+u(j,i)))-(v(j,i-1)*(u(j-1,i-1)+u(j,i-1))))/(4*dx));
        b3=(v(j+1,i)-2*v(j,i)+v(j-1,i))/(dy^2);            %solving N-S for v Velocity
        b4=(v(j,i+1)-2*v(j,i)+v(j,i-1))/(dx^2);
          
        B=b1+(b3+b4)/ratio;
      
        v1(j,i)=v(j,i)+dt*(B-(p(j,i-1)-p(j-1,i-1))/dy);
       
    end
end


%apply boundary conditions 
u1(Nx/2-Nx/5+2:Nx/2+Nx/5,Ny/2-Ny/5+2:Ny/2+Ny/5) = 0.0;                   
v1(Nx/2-Nx/5+2:Nx/2+Nx/5,Ny/2-Ny/5+2:Ny/2+Ny/5) = 0.0;


%iterate for pressure and velocity corrections
 for iteration=1:total_iterations           % Iteration Loop   

for j=1:Ny
    for i=1:Nx
        
        residual1(j,i)=(u1(j,i+1)-u1(j,i)+v1(j+1,i)-v1(j,i))/(-J_a*dt);    %calculate residuals from continuity
   
    end
end


for j=1:Ny
    for i=1:Nx
        
        residual(Nx*(j-1)+i,1)=residual1(j,i);                          %converting residual from a matrix to a vector
    
    end
end

dp=J\residual;                                              %changes in pressure field

for j=1:Ny
    for i=1:Nx
        dp1(j,i)=dp(Nx*(j-1)+i,1);                             %converitng changes in pressure field from a vector to a matrix
    end
end

for j=2:Ny
    for i=2:Nx
        u1(j,i)=u1(j,i)+relaxation_factor*(dp1(j,i-1)-dp1(j,i))*dt/dx;                 %u velocity correction
    end
end

for j=2:Ny
    for i=2:Nx+1
        v1(j,i)=v1(j,i)+relaxation_factor*(dp1(j-1,i-1)-dp1(j,i-1))*dt/dy;             %v velocity correction
    end
end


p = p + relaxation_factor*dp1;                                      %pressure field correction
u = u1;
v = v1;

residual_max(iteration) = max(abs(residual));                  %output maximum value of residual

if residual_max(iteration) < 1.0e-4                            %stop on convergance
    break
end

 end                  %iterations ends


figure(1)
contourf (p)  %plot pressure field
hold on

UU = u(2:Ny-1,3:Nx);  %select u velocity field (adjust for staggered grid)
VV = v(2:Ny-1,3:Nx);  %select v velocity field (adjust for staggered grid)

[X,Y]=meshgrid(2:1:Nx-1,2:1:Ny-1);   %vector plot
q=quiver(X,Y,UU,VV,1);
q.Color = 'black';
axis equal;

%draw a square
v_ = [Nx/2-Nx/5+1 Ny/2-Ny/5+1.5;  Nx/2-Nx/5+1 Ny/2+Ny/5+0.5; Nx/2+Nx/5 Ny/2-Ny/5+1.5; Nx/2+Nx/5 Ny/2+Ny/5+.5];
f_ = [2 1 3 4];
p_=patch('Faces',f_,'Vertices',v_,'FaceColor','red');
p_.EdgeColor='none';
p_.FaceColor='white';
xlabel ('x-dimension')
ylabel ('y-dimension')
title ('Pressure and velocity field around the square')
colorbar

figure()
plot(residual_max)  %residual plot
xlabel ('Iteration number')
ylabel ('Maximum value of residual')
title ('Convergence plot')



