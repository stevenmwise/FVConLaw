function [cM,cP] = WENO5(c)
%
% Fifth-Order WENO Reconstructions:
%
mx = length(c);
cM = zeros(mx,1);
cP = zeros(mx,1);
%
a11 =  02.0/6.0; a12 = -07.0/6.0; a13 =  11.0/6.0;
a21 = -01.0/6.0; a22 =  05.0/6.0; a23 =  02.0/6.0;
a31 =  02.0/6.0; a32 =  05.0/6.0; a33 = -01.0/6.0;
g1 = 01.0/10.0;
g2 = 06.0/10.0;
g3 = 03.0/10.0;
epsln = 1.0e-06;
%
% Periodic boundary conditions:
cim2 = circshift(c, 2);
cim1 = circshift(c, 1);
cip1 = circshift(c,-1);
cip2 = circshift(c,-2);
cip3 = circshift(c,-3);
%
% 3 quadratic polynomials:
f1(1:mx) = a11*cim2(1:mx)+a12*cim1(1:mx)+a13*c   (1:mx);
f2(1:mx) = a21*cim1(1:mx)+a22*c   (1:mx)+a23*cip1(1:mx);
f3(1:mx) = a31*c   (1:mx)+a32*cip1(1:mx)+a33*cip2(1:mx);
%
% Smoothness measures:
b1(1:mx) ...
  = 13/12*(    cim2(1:mx)-2.0*cim1(1:mx)+    c   (1:mx)).^2 ...
  + 01/04*(    cim2(1:mx)-4.0*cim1(1:mx)+3.0*c   (1:mx)).^2;
b2(1:mx) ...
  = 13/12*(    cim1(1:mx)-2.0*c   (1:mx)+    cip1(1:mx)).^2 ...
  + 01/04*(    cim1(1:mx)               -    cip1(1:mx)).^2;
b3(1:mx) ...
  = 13/12*(    c   (1:mx)-2.0*cip1(1:mx)+    cip2(1:mx)).^2 ...
  + 01/04*(3.0*c   (1:mx)-4.0*cip1(1:mx)+    cip2(1:mx)).^2;
%
% Weights:
v1(1:mx) = g1./(epsln+b1(1:mx)).^2;
v2(1:mx) = g2./(epsln+b2(1:mx)).^2;
v3(1:mx) = g2./(epsln+b3(1:mx)).^2;
w1(1:mx) = v1(1:mx)./(v1(1:mx)+v2(1:mx)+v3(1:mx));
w2(1:mx) = v2(1:mx)./(v1(1:mx)+v2(1:mx)+v3(1:mx));
w3(1:mx) = v3(1:mx)./(v1(1:mx)+v2(1:mx)+v3(1:mx));
%
% 5th-order WENO reconstruction c_{i+1/2}^-
cM(1:mx) = f1(1:mx).*w1(1:mx) ...
         + f2(1:mx).*w2(1:mx) ...
         + f3(1:mx).*w3(1:mx);
%
% 3 quadratic polynomials:
f1(1:mx) = a11*cip3(1:mx)+a12*cip2(1:mx)+a13*cip1(1:mx);
f2(1:mx) = a21*cip2(1:mx)+a22*cip1(1:mx)+a23*c   (1:mx);
f3(1:mx) = a31*cip1(1:mx)+a32*c   (1:mx)+a33*cim1(1:mx);
%
% Smoothness measures:
b1(1:mx) ...
  = 13/12*(    cip3(1:mx)-2.0*cip2(1:mx)+    cip1(1:mx)).^2 ...
  + 01/04*(    cip3(1:mx)-4.0*cip2(1:mx)+3.0*cip1(1:mx)).^2;
b2(1:mx) ...
  = 13/12*(    cip2(1:mx)-2.0*cip1(1:mx)+    c   (1:mx)).^2 ...
  + 01/04*(    cip2(1:mx)               -    c   (1:mx)).^2;
b3(1:mx) ...
  = 13/12*(    cip1(1:mx)-2.0*c   (1:mx)+    cim1(1:mx)).^2 ...
  + 01/04*(3.0*cip1(1:mx)-4.0*c   (1:mx)+    cim1(1:mx)).^2;
%
% Weights:
v1(1:mx) = g1./(epsln+b1(1:mx)).^2;
v2(1:mx) = g2./(epsln+b2(1:mx)).^2;
v3(1:mx) = g2./(epsln+b3(1:mx)).^2;
w1(1:mx) = v1(1:mx)./(v1(1:mx)+v2(1:mx)+v3(1:mx));
w2(1:mx) = v2(1:mx)./(v1(1:mx)+v2(1:mx)+v3(1:mx));
w3(1:mx) = v3(1:mx)./(v1(1:mx)+v2(1:mx)+v3(1:mx));
%
% 5th-order WENO reconstruction c_{i+1/2}^+
cP(1:mx) = f1(1:mx).*w1(1:mx) ...
         + f2(1:mx).*w2(1:mx) ...
         + f3(1:mx).*w3(1:mx);
end