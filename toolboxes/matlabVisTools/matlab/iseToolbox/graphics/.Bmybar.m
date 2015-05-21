function mybar(f,s,xstr,ystr,col,LineWidth,w,we);
%
% mybar(f,[s],[xstr],[ystr],[color],[linewidth],[barwidth],[errorbarwidth]);
% 
%  Draws a bar graph with heights determined by f
%  and error-bars determined by s.
%
%  if input matrices are 2d, the ROWS are plotted as
%  separate sub-bars in different colors.
%
%  xlabels set to 'xstr' (columns)
%  legend set to 'ystr' (rows)
%
%  color is a string of colors, one for each subbar
%  (default is 'brgmyc')
%
%  linewidth (default 2)
%
%  barwidth (default 0.65)
%  errorbarWidth (default 0.3)
%
%  Example:
%    f = rand(5,1)+1;
%    s = rand(5,1)/2;
%    xstr = str2mat('one','two','three','four','five');
%    ystr = '';
%    mybar(f,s,xstr,ystr,'g');

%  7/11/95  gmb Wrote it.
%  4/11/96  gmb Revised it to plot 2d data sets
%  11/19/97 gmb Converted it to Matlab version 5.0
%  12/30/97 djh Fixed bug in xlabels, for some reason it was
%               putting the 1st label under the 2nd bar, and so on.
%  1/24/97  djh Added linewidth as optional input arg

if ~exist('col')
  col='brgmyc';
end
if ~exist('LineWidth')
  LineWidth=2;
end
if ~exist('w')
  w=0.65;
end
if ~exist('we')
  we=0.3;
end

if nargin==1
  s=zeros(size(f));
end

if(size(f,1)==1)
  f=f';
  s=conj(s');
end

if (isreal(s))
  s = s+sqrt(-1)*s;
end

nsubbars=size(f,2);
we=we/nsubbars;

heldstate=ishold;
plot(0,0)
set(gca,'XLim',[1-w,size(f,1)+w]);
    
yhilim = max(max(f+s))*1.2;
if ~isnan(yhilim)
  %set(gca,'YLim',[0,yhilim]);
end

subx=linspace(-w/2,w/2,nsubbars+1);

for subbar=1:nsubbars
  x=[];
  y=[];
  for i=1:size(f,1);
    x=[x';[i+subx(subbar),i+subx(subbar),i+subx(subbar+1),i+subx(subbar+1)]]';
    y=[y';[0,f(i,subbar),f(i,subbar),0]]';
  end
  
  fill(x,y,col(subbar));
  
  a=line(get(gca,'XLim'),[0,0]);
  set(a,'Color','k');
	
  errx=mean([subx(subbar),subx(subbar+1)]);
  for i=1:size(f,1)
    line([i+subx(subbar),i+subx(subbar),i+subx(subbar+1),i+subx(subbar+1)], ...
	[0,f(i,subbar),f(i,subbar),0],'Color','k','LineWidth',LineWidth);
    g=line([i+errx,i+errx],[f(i,subbar)-imag(s(i,subbar)),f(i,subbar)+real(s(i,subbar))],...
	'Color','k','LineWidth',LineWidth);
    g=line([i+errx-we/2,i+errx+we/2],[f(i,subbar)-imag(s(i,subbar)),f(i,subbar)-imag(s(i,subbar))],...
	'Color','k','LineWidth',LineWidth);
    g=line([i+errx-we/2,i+errx+we/2],[f(i,subbar)+real(s(i,subbar)),f(i,subbar)+real(s(i,subbar))],...
	'Color','k','LineWidth',LineWidth);
  end
  hold on
end 					%subbars

% Xlabels
set(gca,'xLimMode','manual');
set(gca,'XTick',[1:size(f,1)]);
if (exist('xstr'))
  if ~isempty(xstr)
    set(gca,'XTickLabel',xstr);
  end
end

% Legend
if (exist('ystr'))
  if ~isempty(ystr)
    %heighten the graph
    ylim=get(gca,'Ylim');
    set(gca,'YLim',[ylim(1),ylim(2)*1.5]);
    bs=8;
    
    posgca=get(gca,'Position');
    posgcf=get(gcf,'Position');	
    dx=diff(get(gca,'XLim'))/(posgca(3)*posgcf(3));
    dy=diff(get(gca,'YLim'))/(posgca(4)*posgcf(4));
    
    bw=bs*dx; 				%box width
    bh=bs*dy; 				%box height	
    
    yspace=3.5*bh; 			%vertical spacing	
    
    xlim=get(gca,'Xlim');
    ylim=get(gca,'Ylim');
    
    xc=xlim(1)+30*dx;
    yc=ylim(2)-30*dy;
    bx=[xc-bw,xc+bw,xc+bw,xc-bw,xc-bw];
    by=[yc+bh,yc+bh,yc-bh,yc-bh,yc+bh];
    
    textx=xc+bw+10*dx;
    texty=yc;
    
    %draw box and text
    for i=1:size(ystr,1);
      fill(bx,by,col(i));
      line(bx,by,'LineWidth',LineWidth,'Color','k');
      text(textx,texty,ystr(i,:));
      by=by-yspace;
      texty=texty-yspace;
      
    end
  end
end
	
if (heldstate==0)
	hold off
end


