function [u] = SSPRK3(uo,dx,dt,sigma,reconFunHndl,numFluxHndl)
%
% first RK stage
[uM,uP] = reconFunHndl(uo);
fR = numFluxHndl(uM,uP,dx,dt,sigma);    
fL = circshift(fR,1);
u = uo-dt/dx*(fR-fL);
%
% second RK stage
[uM,uP] = reconFunHndl(u);
fR = numFluxHndl(uM,uP,dx,dt,sigma);    
fL = circshift(fR,1);
u = 3/4*uo+1/4*(u-dt/dx*(fR-fL));
%
% third RK stage
[uM,uP] = reconFunHndl(u);
fR = numFluxHndl(uM,uP,dx,dt,sigma);    
fL = circshift(fR,1);
u = 1/3*uo+2/3*(u-dt/dx*(fR-fL));
%
end