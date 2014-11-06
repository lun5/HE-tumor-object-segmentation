function myplot(x,y,str,sz,LineWidth)

%MYPLOT	myplot(x,y,str,[sz],[LineWidth])
%
%Plot vectors or matrices with an extended symbol set
%	MYPLOT(X,Y) plots vector X versus vector Y. If X or Y is a matrix,
%	then the vector is plotted versus the rows or columns of the matrix,
%	whichever line up. 
%
%	MYPLOT(Y) plots the columns of Y versus their index.
%	If Y is complex, MYPLOT(Y) is equivalent to MYPLOT(real(Y),imag(Y)).
%	In all other uses of MYPLOT, the imaginary part is ignored.
%
%	Various line types, plot symbols and colors may be obtained with
%	PLOT(X,Y,S) where S is a 1, 2 or 3 character string made from
%	the following characters:
%
%	       y     yellow        .     point		d	diamond
%	       m     magenta       +     plus		s 	square
%	       c     cyan          x     x-mark		t	triangle
%	       r     red           			o	my circle
%	       g     green         :  	 dotted		#	star
%	       b     blue          *     star		a	cycle through (all) symbols
%	       w     white         -     solid		f	fill any of the above
%	       k     black         -.    dashdot	n	not fill
% 	                           --    dashed
%                             	 
%	EXAMPlES:	
%	MYPLOT(X,Y,'c+') plots a cyan plus at each data point.
%	MYPLOT(x,y,'cdf-') plots filled cyan diamonds with interpolating lines.
%	MYPLOT(x,y,'fs--') plots filled squares with interpolating dashed lines.
%		colors are cycled for each column of y.
%	MYPLOT(x,y,'a') cycles through symbols (first filled 'stdc#', then not filled)
%	MYPLOT(x,y,'an') cycles through open symbols
%	MYPLOT(x,y,'awf-') cycles through white filled symbols and interpolates solid lines.
%
%	The MYPLOT command, if no color is specified, makes automatic use of
%	the colors specified by the axes ColorOrder property.  The default
%	ColorOrder is listed in the table above for color systems where the
%	default is yellow for one line, and for multiple lines, to cycle
%	through the first six colors in the table.  For monochrome systems,
%	MYPLOT cycles over the axes LineStyleOrder property.
%
%	MYPLOT returns a column vector of handles to LINE objects, one
%	handle per line. (handles to new symbols are not returned)
%
%	See also SEMILOGX, SEMILOGY, LOGLOG, GRID, CLF, CLC, TITLE,
%	XLABEL, YLABEL, AXIS, AXES, HOLD, PLOT, and SUBPLOT MYERRORBAR MYBAR.

%	written 9/27/95	gmb

if nargin<3
	str='-';
end

if ~exist('LineWidth');
  LineWidth=1;
end

symlist='stdo#a';

% parameter values for symbols
sa=[45,0,0,0,0]; 			%starting angle
da=[90,120,90,11.25,144];               %delta-angle
na=[4,3,4,32,5];                        %number of corners
sc=[1,1,1,0.75,1];                      %scale factor
if (size(y,1)==1)
	y=y';
end
collist='ymcrgbwk';                     

if (nargin<=3)
	sz=6;
end

hld=get(gca,'NextPlot');

%look for special characters
sym=0;
for i=1:length(symlist);
	k=(find(symlist(i)==str));
	if k>0
		sym=i;
		str=[str(1:k-1)';str(k+1:length(str))']';
		i=i-1;
	end
end
%fill region?
id=find(str=='f');
ptch=0;
if id>0
	str=[str(1:id-1)';str(id+1:length(str))']';
	p=1;
	ptch=1;
end
id=find(str=='n');
if id>0
	str=[str(1:id-1)';str(id+1:length(str))']';
	p=0;
	ptch=1;
end
%fix annoying way MATLAB skips values with NAN's
ynonan=y;
nanstofix=isnan(y);
if (sum(sum(nanstofix))>0)
	nanstofix(1,:)=zeros(size(y(1,:)));
	nanstofix(size(y,1),:)=zeros(size(y(1,:)));
	nanstofix=find(nanstofix);
	for i=nanstofix'
		x0=x(mod(i,length(x))-1);
		x1=x(mod(i,length(x)));
		x2=x(mod(i,length(x))+1);
		y0=y(i-1);
		y2=y(i+1);
		ynonan(i)=y0+(y2-y0)*(x1-x0)/(x2-x0);
	end
end

if sym==0
	plot(x,ynonan);
else
	%determine color
	col='y';
	k=0;
	for i=1:length(str);
		k=find(str(i)==collist);
		if k>0
			col=collist(k);
		end
	end	
	%determine color
	cyccolor=1;
	for i=1:length(str);
		k=(find(str(i)==collist));
		if k>0
			col=collist(k);
			cyccolor=0;
		end
	end	
	%reset axes if x or y are outside of XLim or YLim
	xlim=get(gca,'XLim');
	ylim=get(gca,'YLim');
	


	%plot non-special portion of graph
	if (length(str)>0 & ~(length(str)==1 & ~cyccolor));
		plot(x,ynonan,str);
		hold on
	end
	set(gca,'XLim',[min(mmin(x(~isnan(x))),xlim(1)),max(mmax(x(~isnan(x))),xlim(2))]);
	set(gca,'YLim',[min(mmin(y(~isnan(y))),ylim(1)),max(mmax(y(~isnan(y))),ylim(2))]);


	%cycle symbols?
	cycsymb=0;
	if (sym==(find(symlist=='a'))) 
		cycsymb=1;
	end
	%determine symbol size
	posgca=get(gca,'Position');
	posgcf=get(gcf,'Position');	
	dx=sz*diff(get(gca,'XLim'))/(posgca(3)*posgcf(3));
	dy=sz*diff(get(gca,'YLim'))/(posgca(4)*posgcf(4));
 	for j=1:size(y,2)
		if (cyccolor==1)
			col=collist(mod(j-1,length(collist)-2)+1);
		end
		if (cycsymb==1)
			sym=mod(j-1,length(symlist)-1)+1;
			if(ptch==0)
				p=floor(mod(j-1,(length(symlist)-1)*2)/(length(symlist)-1));
			end
		end
		xlist=[];
		ylist=[];
		for i=1:size(y,1)
			a=[sa(sym)*pi/180:da(sym)*pi/180:sa(sym)*pi/180+(na(sym)+1)*da(sym)*pi/180];
			xlist=[xlist;x(i)+sin(a)*dx*sc(sym)];
			ylist=[ylist;y(i,j)+cos(a)*dy*sc(sym)];
		end

		if (ptch)
			patch(xlist',ylist',col,'Clipping','off');
		end
		line(xlist',ylist','Color',col,'Clipping','off','LineWidth',LineWidth);	
	end
end
if(hld(1)=='r')
	hold off
end


