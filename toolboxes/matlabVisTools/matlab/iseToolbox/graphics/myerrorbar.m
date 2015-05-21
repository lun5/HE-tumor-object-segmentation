function myerrorbar(x,y,e,col);
% MYERRORBAR
%	myerrorbar(x,y,e,[c]);
%	Draws errorbars on the current plot. 
%	Works best following 'plot(x,y)'
%	Also Good for overlaying data points on model fit.
%	myerrorbar differs from 'errorbar' in two ways:
%	1. User can define color of errorbar with 'c';
%	   By default, colors of errorbars are same as current plot.
%	2. Lines between points are not drawn.
%
% Example:
%	x=linspace(0,1,10);
%	m=x.^2;
%	sd=rand(size(x))/5;
%	myplot(x,m,'fo-')
%	myerrorbar(x,m,sd) 

%	4/16/95 gmb	Wrote it.


if (size(x,1)==1)
	x=x';
end
if (size(y,1)==1)
	y=y';
end
if (size(e,1)==1 | size(e,1) ==2)
	e=e';
end

if isreal(e)
  e = e+sqrt(-1)*e;
end

if (size(x,2)==1 & size(y,2)>1)
	x=x*ones(1,size(y,2));
end

c=get(gca,'ColorOrder');
hld=get(gca,'NextPlot');

xs=1; %width of errorbar (percent)

xlim=get(gca,'XLim');
dx=(xlim(2)-xlim(1))*xs/100;
hold on
for i=1:size(x,2);
	if (nargin<=3)
		col=c(i,:);
	end
	for j=1:size(x,1);
	   if (e(j,i)>0)

		h=line([x(j,i)-dx,x(j,i)+dx],[y(j,i)-real(e(j,i)),y(j,i)-real(e(j,i))],'Color',col);
		h=line([x(j,i)-dx,x(j,i)+dx],[y(j,i)+imag(e(j,i)),y(j,i)+imag(e(j,i))],'Color',col);
		h=line([x(j,i),x(j,i)],[y(j,i)-real(e(j,i)),y(j,i)+imag(e(j,i))],'Color',col);
	   end
       end
end
if(hld(1)=='r')
	hold off
end

