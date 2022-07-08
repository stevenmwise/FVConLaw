function u = burgers(b,xo,t)
%
mx = length(xo);
x = zeros(mx,1);
%
for j = 1:mx
  x(j) = fzero(@(s)myFunc(s,xo(j),t,b),0);
end
%
u = b*sin(x);
%
end