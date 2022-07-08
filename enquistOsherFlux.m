function fHat = enquistOsherFlux(uM,uP,dx,dt,sigma)
%
mx = length(uM);
fHat = zeros(mx,1);
%
for i=1:mx
  if (uM(i)>0)
    a = fluxFunction(uM(i));
  else
    a = 0;
  end
  if (uP(i)>0)
    b = 0;
  else
    b = fluxFunction(uP(i));
  end
  fHat(i) = a+b;
end
%
return