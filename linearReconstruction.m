function [uM,uP] = linearReconstruction(u)
%
mx = length(u);
uM = zeros(mx,1);
uP = zeros(mx,1);
%
um1 = circshift(u,+1);
up1 = circshift(u,-1);
sigma = minMod2(u-um1,up1-u);
%
% The value of u on the left (minus) side of i+1/2:
uM = u+sigma/2.0;
%
% The value of u on the right (plus) side of i+1/2:
uP = u-sigma/2.0;
uP = circshift(uP,-1);
%
end