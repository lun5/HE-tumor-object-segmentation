function widen(p)
%widen([p]) 
%spreads x-axis limits by a proportion p>1. Designed to pull 
%extreme data points off of the vertical axes. 
%default: p=1.1;

%2/27/97 gmb	wrote the three lines of code.

if (nargin==0)
	p=1.1;
end

XLim=get(gca,'XLim');
XLim=([-.5,0.5]*diff(XLim)*p)+mean(XLim);
set(gca,'XLim',XLim);
