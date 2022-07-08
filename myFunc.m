function [f,fp] = myFunc(xo,x,t,b)

f  = x - xo - (b*sin(xo)).*t;
fp = -1 - t.*b.*cos(xo);

return