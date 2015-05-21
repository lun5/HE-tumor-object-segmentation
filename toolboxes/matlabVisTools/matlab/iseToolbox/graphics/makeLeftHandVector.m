function p = makeLeftHandVector(v,p)
%
%AUTHOR:  Engel, Wandell
%DATE:    03.22.95
%
% This routine takes two PERPENDICULAR vectors, v, and p.  It decides
% whether p is 90 deg counter-clockwise to v.  If it is not, but rather it is
% 90 deg clockwise, then it returns a vector that is flipped by 180 deg.
%

if v'*p> 3*eps
  error('Input vectors must be perpendicular')
end

aV = atan2(v(2),v(1));
aP = atan2(p(2),p(1));
if aV < 0
 aV = 2*pi + aV;
end
if aP < 0
 aP = 2*pi + aP;
end

diffAngle = aP - aV;
if abs(diffAngle) - pi/2 < 3*eps;
  s = sign(diffAngle);
elseif diffAngle - (-6*pi/4) < 3*eps
  s = 1;
elseif diffAngle - (6*pi/4) < 3*eps
  s = -1;
else
 error('makeLeftHandVector:  Routine bug error')
end
p = s*p;

