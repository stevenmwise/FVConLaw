function [vv] = exactSolution(t,x,alpha,beta)

% location of shock
q = pi+alpha*t;
if (q > 2*pi)
  d = q-2*pi*floor(q/(2*pi));
  x = x+2*pi*floor(q/(2*pi));
else
  d = q;
end

d2 = q;

x1 = x(x< d2);
x2 = x(x>=d2);

if (~isempty(x1))
%
% coordinate transformation
  xn = x1-alpha*t;
%
% solve burgers using the coordinate transformation
  u = burgers(beta,xn,t);
  v = alpha+u;
%
  x1 = x1-2*pi*floor(q/(2*pi));
    
  if (~isempty(x2))
    x2 = flip(2*d2-x2);
    xn2 = x2-alpha*t;

    u2 = burgers(beta,xn2,t);
    v2 = alpha+u2;

    x2 = x2-2*pi*floor(q/(2*pi));
    xx = [x1; flip(2*d-x2)];

    vv = [v; flip(2*alpha-v2)];
  else
    xx = x1;
    vv = v;
  end
elseif (~isempty(x2))
  x2 = flip(2*d2-x2);
  xn2 = x2-alpha*t;

  u2 = burgers(beta,xn2,t);
  v2 = alpha+u2;

  x2 = x2-2*pi*floor(q/(2*pi));
  xx = flip(2*d-x2);
  vv = flip(2*alpha-v2);
end

end


