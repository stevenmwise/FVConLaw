function fHat = laxWendFlux(uM,uP,dx,dt,sigma)
%
mx = length(uM);
fHat = zeros(mx,1);
%
fM = fluxFunction(uM); 
fP = fluxFunction(uP);
%
fHat = 0.5*(fM+fP-(dt/dx)*0.5*(uM+uP).*(fP-fM));
%
return