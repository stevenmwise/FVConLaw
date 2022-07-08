function fHat = roeFlux(uM,uP,dx,dt,sigma)
%
fM = fluxFunction(uM);
fP = fluxFunction(uP);
%
mx = length(uM);
s    = zeros(mx,1);
fHat = zeros(mx,1);
%
for i = 1:mx
  if (uP(i) == uM(i))
    s(i) = uM(i);
  else
    s(i) = (fP(i)-fM(i))/(uP(i)-uM(i));
  end
end
%
fHat = 0.5*(fM+fP-sign(s).*(fP-fM));
%
return