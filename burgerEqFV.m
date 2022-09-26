%
% This script is used to obtain approximations to the solution of
% inviscid Burger's Equation:
%
% u_t+(f(u))_x = 0
%
% Periodic boundary conditions are assumed over the domain [a,b].
%
% We use a finite volume discretization method with either a
% 5th-order WENO reconstruction or a linear slope-limited
% reconstruction of the function. Numerical fluxes can be chosen
% via the flag numerFluxNum from the following options:
%
% numerFluxNum = 1: Roe
% numerFluxNum = 2: Godunov
% numerFluxNum = 3: LaxFriedrichs
% numerFluxNum = 4: Enquist-Osher
% numerFluxNum = 5: Lax Wendroff
%
% Reconstructions can be chosen via the flag reconFunNum:
%
% reconFunNum = 1: WENO5
% reconFunNum = 2: Linear reconstruction with minmod limiter
%
% Time integration is carried out by strongly stability preserving
% (SSP) Runge Kutta methods. The scheme can be chosen via the flag 
% timeIntNum from the following list of options:
%
% timeIntNum = 2: second-order SSP Runge-Kutta method
% timeIntNum = 3: third-order SSP Runge-Kutta method
%
clear all; close all; clc;
%
% Flux selection:
numerFluxNum = 4;
numerFluxHndl = numerFluxSelection(numerFluxNum);
%
% Reconstruction selection:
reconFunNum = 1;
reconFunHndl = reconFunSelection(reconFunNum);
%
% Time integration selection:
timeIntNum = 2;
timeIntHndl = timeIntSelection(timeIntNum);
%
% Spatial domain:
a = 0;
b = 2*pi;
%
% Initial condition parameters: u(x,0) = alpha + beta*sin(x)
alpha = 0.0;
beta  = 1.0;
%
% Number of finite volume cells:
mx  = 100;
%
% Final time:
finalTime = 1.2;
% 
% Edge grid points:
x = linspace(a,b,mx+1)';  
%
% Grid spacing:
dx = (b-a)/(mx);   
%
% Initial data:
uo = initialData(alpha,beta,mx,dx,x);
%
% Initial time step:
dt = dx/(2*max(max(uo)));
%
currentTime = 0.0;
%
% Time integration loop:
while (currentTime < finalTime)
  sigma  = max(abs(uo));
%
  u = timeIntHndl(uo,dx,dt,sigma,reconFunHndl,numerFluxHndl);
%
  uo = u;
  currentTime = currentTime+dt;
%
  dt = min(finalTime-currentTime,dx/(2*max(max(uo))));
end
%
currentTime
%
mxFine = 300;
xFine = linspace(a,b,mxFine+1)'; 
uFine = zeros(mxFine,1);
%
dxFine = (b-a)/mxFine;
%
for i=1:mxFine
  uFine(i) = 1/dxFine*integral(@(s)exactSolution(currentTime,s', ...
    alpha,beta),xFine(i),xFine(i+1),'ArrayValued',true);
%  uFine(i) = 1/dxFine*integral(@(s)exactSolution(currentTime,s', ...
%    alpha,beta),xFine(i),xFine(i+1));
end  
%
% plot results
figure
xCntr(1:mx) = x(1:mx)+dx/2;
plot(xCntr,u,'o')
hold on;
xFineCntr(1:mxFine) = xFine(1:mxFine)+dxFine/2;
plot(xFineCntr,uFine)
hold off;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Embedded functions below:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function numerFluxHndl = numerFluxSelection(numerFluxNum)
%
switch numerFluxNum
  case 1
    numerFluxHndl = @roeFlux;
  case 2
    numerFluxHndl = @godunovFlux;
  case 3
    numerFluxHndl = @laxFriedFlux;
  case 4
    numerFluxHndl = @enquistOsherFlux;
  case 5
    numerFluxHndl = @laxWendFlux;
end
end
%
function reconFunHndl = reconFunSelection(reconFunNum)
%
switch reconFunNum
  case 1
    reconFunHndl = @WENO5;
  case 2
    reconFunHndl = @linearReconstruction;
end
end
%
function timeIntHndl = timeIntSelection(timeIntNum)
%
switch timeIntNum
  case 2
    timeIntHndl = @SSPRK2;
  case 3
    timeIntHndl = @SSPRK3;
end
end
%
function uo = initialData(alpha,beta,mx,dx,x)
%
% Cell centered intial data:
uo = zeros(mx,1); 
%
% Compute cell averages at t = 0:
for i=1:mx
%  uo(i) = 1/dx*integral(@(s)alpha+beta*sin(s),x(i),x(i+1), ...
%    'ArrayValued',true);
  uo(i) = 1/dx*integral(@(s)alpha+beta*sin(s),x(i),x(i+1));
end
end
