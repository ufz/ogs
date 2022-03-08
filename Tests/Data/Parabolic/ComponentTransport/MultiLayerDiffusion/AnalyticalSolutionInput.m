% This Matlab script is used to produce the semi-analytical solution of the two-layer diffusion problem. If you want to reproduce this semi-analytical solution locally, you need to download the Matlab code from the MultDiff project on github (see https://github.com/elliotcarr/MultDiff).

close all, clear all, clc

% Add directory to current path
addpath('..')

% Parameters
m     = 3;                                                        % Number of layers
kappa = [5.555555555555556e-10, 8.333333333333333e-11 8.333333333333333e-11];           % Diffusivities
l0    = 0.0;                                                      % Left end of slab
lm    = 20.0;                                                     % Right end of slab
l     = [0.625, 15];                                                  % Location of interfaces
u0    = @(x) zeros(size(x));                                      % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,1.0};                                % Boundary condition (x = 0 m)
Rbnd  = {'Neumann',0.0,1.0,0.0};                                  % Boundary condition (x = 20 m)
tspan = [31536000000,315360000000,3153600000000,31536000000000];  % Times at which to compute solution

%% Compute analytical solution
options.N = 500; % Number of eigenvalues
options.NX = 250; % Number of divisions in each slab
%options.prob = 'C';

[u,x] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect', options);

% Plot
figure;
for i = 1:m-1
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
M=[u x];
writematrix(M, 'M.csv');
plot(x,u,'b','LineWidth',2.0);
axis([0,20,-0,1.1]);
xlabel('$x$','Interpreter','LaTeX','FontSize',20);
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20);
set(gca,'FontSize',14,'Layer','top');
