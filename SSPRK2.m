function [u] = SSPRK2(uo,dx,dt,sigma,reconFunHndl,numFluxHndl)
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
u = 1/2*uo+1/2*(u-dt/dx*(fR-fL));
%
end