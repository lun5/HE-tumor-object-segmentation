function PlotFitDesign(m,sd,x,pred,predx,xlab,ylab,legstr,xlogflag,symbols);
%PlotFitDesign(m,sd,x,pred,predx,[xlab],[ylab],[legstr],[xlogflag],[symbols]);
%
%Plots 'plot(x,m)' with errorbars of size 'sd'.
%Labels axes with 'xlab' and 'ylab'.
%Adds legend with 'legstr'.
%if xlogflag=='y', plots are on log10 x axis.
%
%SEE PlotDesign

%1/16/96 gmb	Wrote it.
%7/28/96 gmb	Created PlotFitDesign from PlotDesign

%plot it

hashflag=0;


if (nargin>8)
	if (xlogflag(1)=='y')
		if(x(1)==0)
			hashflag=1;
			x(1)=1;
		end
		x=log10(x);
		predx=log10(predx);
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
plot(predx,pred,'-');
myerrorbar(x,m,sd)

hold on
if nargin>7
	if (legstr ~= '')
		legend(legstr);
	end
end


set(gca,'XLim',[min(x),max(x)]);
set(gca,'XTick',x);
%ylim=get(gca,'YLim');
%ylim(1)=0;
%set(gca,'YLim',ylim);


if nargin>5
	xlabel(xlab);
	ylabel(ylab);
end
hold off

if (xlogflag(1)=='y')
	logx2raw(10);
end


%Last, but certainly not least...
if nargin>9
  for i=1:size(m,2)
    myplot(x,m(:,i),symbols(i,:))
    hold on
  end
  hold off
else
  myplot(x,m,'af')
end
if(hashflag==1)
  HashMark;
end




