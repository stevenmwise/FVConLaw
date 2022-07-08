function fHat = godunovFlux(uM,uP,dx,dt,sigma)
%
mx = length(uM);
fHat = zeros(mx,1);
%
for i=1:mx
  if (uM(i)*uP(i)<0 && uM(i)<=uP(i))
    fHat(i) = 0;
  elseif (uM(i)*uP(i)<0 && uM(i)>uP(i))
    fHat(i) = max([fluxFunction(uM(i)) fluxFunction(uP(i))]);
  elseif (uM(i)*uP(i)>=0 && uM(i)<=uP(i))
    fHat(i) = min([fluxFunction(uM(i)) fluxFunction(uP(i))]);
  elseif (uM(i)*uP(i)>=0 && uM(i)>uP(i))
    fHat(i) = max([fluxFunction(uM(i)) fluxFunction(uP(i))]);
  end
end

return