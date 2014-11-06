function PlotDesign(m,sd,x,xlab,ylab,legstr,xlogflag,symbols);
%PlotDesign(m,sd,x,[xlab],[ylab],[legstr],[xlogflag],[symbols]);
%
%Plots 'plot(x,m)' with errorbars of size 'sd'.
%Labels axes with 'xlab' and 'ylab'.
%Adds legend with 'legstr'.
%if xlogflag=='y', plots are on log10 x axis.
%
%Example

%x=linspace(0,1,11)';
%m=[x.^2 x.^3];
%sd=rand(size(m))/10;
%PlotDesign(m,sd,x,'X axis','Y axis',str2mat('squared','cubed'),'y');


%1/16/96 gmb	Wrote it.


%plot it
hashflag=0;

if (nargin>6)
	if (xlogflag(1)=='y')
		if(x(1)==0)
			hashflag=1;
			x(1)=1;
		end
		x=log10(x);
		if (hashflag==1)
			if(length(x)>2)
				x(1)=2*x(2)-x(3);
			else
				x(1)=-inf;
			end
		end
	end
else
	xlogflag='n';
end

myplot(x,m)
myerrorbar(x,m,sd)

set(gca,'XLim',[min(x),max(x)]);
set(gca,'XTick',x);
ylim=get(gca,'YLim');
ylim(1)=0;
set(gca,'YLim',ylim);

if nargin>7
  for i=1:size(m,2)
    myplot(x,m(:,i),symbols(i,:))
    hold on
  end
  hold off
else
  myplot(x,m,'af')
end

hold on
if nargin>5	
	if (legstr ~= '')
		legend(legstr);
	end
end

if nargin>3
	xlabel(xlab);
	ylabel(ylab);
end
hold off

if (xlogflag(1)=='y')
	logx2raw(10);
end

if(hashflag==1)
	HashMark;
end


