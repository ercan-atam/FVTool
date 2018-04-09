
%% This is a 1D transient convection-diffusion example with piece-wise heat generation, which is Example 8.3 
% from Versteg_Malalasekera (Computational Fluid Dynamics, 2nd edt, page 258), which also can be found from the following link:
% http://www.skb.se/darcytools/wp-content/uploads/2017/05/DarcyTools_EC_20170427_ex2.pdf
% Note: the analytical formula is taken from the link.

%% Clearing, closing

clc; clear;
close all;

%% thermodynamic and other parameters

rho=1; % kg/m^3
alfa=rho;

k=0.03; % kg/ms

x1=0.6; % m
x2=0.2;% m

a=-200;
b=100;


%% Mesh generation

L=1.5; % m

Nx=5000; % number of cells

m = createMesh1D(Nx, L);

figure(1);
plot(m.cellcenters.x, ones(size(m.cellcenters.x)), 'or', ...
     m.facecenters.x, ones(size(m.facecenters.x)), '-+b');
legend('cell centers', 'face centers');
title('Visualization of  the 1D discretized domain');

%% define the boundary conditions

BC = createBC(m); % creates a boundary condition structure

% left boundary: zero temperature
BC.left.a = 0;
BC.left.b = 1;
BC.left.c = 0;

% right boundary: insulated
BC.right.a = 1;
BC.right.b = 0;
BC.right.c = 0;


%% define the diffusion term

D = createCellVariable(m, k); % assign a constant value of k to diffusivity on each cell
D_face = harmonicMean(D); % average diffusivity value on the cell faces

%% define the convection term

u_val= 2; % m/s
u_face = createFaceVariable(m, u_val);


%%  define the source term

f = createCellVariable(m, 0.0);

dx=L/Nx;
ctr_1=0;
for s=min(m.cellcenters.x):dx:max(m.cellcenters.x) 
    ctr_1=ctr_1+1;
    if (s<x1)
        f.value(ctr_1) = a*s+b;
    elseif (s<x1+x2)
        f.value(ctr_1) = 100*s-80;
    else
        f.value(ctr_1) = 0;
    end
end

%% define initial values

T0_val = 0;
T0 = createCellVariable(m, T0_val, BC); % initial values
T = T0; % assign the T0 value of the cells to the current values


%% loop for transient solution

dt = 0.01; % time step in s
final_t = 5; % final time

for t=dt:dt:final_t
    
     [M_trans, RHS_trans] = transientTerm(T, dt, alfa);
    
     M_conv =  convectionUpwindTerm(u_face);
     
    
     M_diff = diffusionTerm(D_face);
    
    [M_bc, RHS_bc] = boundaryCondition(BC);
    
    M = M_conv+M_trans-M_diff+M_bc;
   
    RHS_f=constantSourceTerm(f);
    
    RHS = RHS_trans+RHS_bc+RHS_f;
   
    T = solvePDE(m,M, RHS);
    
end

%% steady-state analytical solution

% analytical solution parameters
Nf=600;
P=rho*u_val/k;
a0=[(x1+x2)*(a*x1+b)+b*x1]/[L];

C3=a0/(2*k*P);

for i=1:Nf
        an(i)=(2*L/(i^2*pi^2))*[((a*x1+a*x2+b)/x2)*cos(i*pi*x1/L)-a-(a*x1+b)/x2*cos(i*pi*(x1+x2)/L)];
end

for i=1:Nf
        C4(i)=an(i)/[k*(P^2+(i*pi/L)^2)];
end

C2_sum=0;
for i=1:Nf
    C2_sum=C2_sum+(-1)^i*C4(i);
end

C2=-exp(-P*L)*(C3/P+C2_sum);


C1_sum=0;
for i=1:Nf
    C1_sum=C1_sum+C4(i);
end
C1=-C2-C1_sum;



ctr=0;
for x_cor=0:dx:L
    
    ctr=ctr+1;
      
    sum=0;
    for i=1:Nf
        sum=sum+C4(i)*[cos(i*pi*x_cor/L)+(L*P)/(i*pi)*sin(i*pi*x_cor/L)];
    end
    
    T_anal(ctr)=C1+C2*exp(P*x_cor)+C3*x_cor+sum;
    
    sum_s=0;
    for i=1:Nf
       sum_s=sum_s+an(i)*cos(i*pi*x_cor/L);
    end
    
    Sr(ctr)=a0/2+sum_s;
    
end     
    
%% comparing results

figure(2)
plot(m.cellcenters.x, T.value(2:end-1));
hold on
plot(m.cellcenters.x, T_anal(1:Nx)','r');
title('Comparing analytical result with the result of FVTool.')
legend('FVTool','Analytical')

