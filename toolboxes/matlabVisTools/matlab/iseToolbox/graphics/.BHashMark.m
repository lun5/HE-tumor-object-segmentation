function HashMark(val,FontSize,FontWeight)
%HASHMARK
%
% HashMark([val],[FontSize],[FontWeight]) 
% places hash marks on x-axis and replaces the
% lower x-limit string with *val*.  -infinity is
% the default.
%
% Example:
%  plot(1:10,1:10)
%  HashMark

% gmb 9/28/95 wrote it.
% gmb 11/7/96 added FontWeight option. 
%              if FontWeight is not given, the graph's axis 
%	      FontWeight is used.

sz=5;
if (nargin==0) 
	val=-inf;
end

if (nargin<2)
	FontSize=get(gca,'FontSize');
end

if (nargin<3)
  FontWeight = get(gca,'FontWeight');
end

posgca=get(gca,'Position');
posgcf=get(gcf,'Position');	
dx=diff(get(gca,'XLim'))/(posgca(3)*posgcf(3));
dy=diff(get(gca,'YLim'))/(posgca(4)*posgcf(4));

xlab=get(gca,'XTickLabels');
x1=min(get(gca,'XLim'));
x=x1+dx*25;
blnk='        ';
xlab(1,:)=blnk(1:length(xlab(1,:)));
set(gca,'XTickLabel',xlab);
ylim=get(gca,'YLim');
y=ylim(1);
for i=x-sz*dx/2:sz*dx:x+sz*dx/2
	a=line([i,i+sz*dx],[y-sz*dy,y+sz*dy],'Clipping','off','Color','k');
end
xc=min(get(gca,'XLim'));
yc=ylim(1)-dy*18;

if (val == -inf)
	text(xc+dx*5*12/FontSize,yc,'8','Rot',90,'FontSize',FontSize,'FontWeight',FontWeight);
	text(xc-dx*10*12/FontSize,yc+dy*6*12/FontSize,'-','FontSize',FontSize,'FontWeight',FontWeight);
elseif finite(val)
	text(xc,yc+dy*6*12/FontSize,num2str(val),'FontSize',FontSize,'FontWeight',FontWeight);
end
set(gca,'YLim',ylim);







