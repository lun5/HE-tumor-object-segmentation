function LabelGraph(var,varname,headername,fontsize)
% LabelGraph(var,varname,[headername],[fontsize])
% Labels the current graph with
%
%	`var(1) = varname(1,:)'
%	`var(2) = varname(2,:)'
%	     .
%	     .
%	     .
%	`var(n) = varname(n,:)'
%
% in the upper-left corner of the graph.
% Any extra strings sent are interpreted as a header or title to the label.
% Any string preceeded by `#' is interpreted as a symbol. (i.e. '#s' is a sigma)
%
% Designed to label graphs of model fits with parameter values.
%
% Example:
%  plot(0:10,1+(0:10)/2)
%  LabelGraph([0.5,1],str2mat('Slope','Y-intercept'))
%
% Can also place text in the upper left: 
% LabelGraph([],[],str2mat('First string','Second string'))

%1/17/96	gmb	wrote it.



if (nargin<2)
	error('Not enough input arguments');
	return
end

if nargin<4
  fontsize=12;
end

k=0.05;

YLim=get(gca,'YLim');
XLim=get(gca,'XLim');

xc=XLim(1)+diff(XLim)*k;
yc=YLim(2)-diff(YLim)*k;
dy=diff(YLim)*k;
spx=diff(XLim)*.02;
%check for header
if nargin>=3;
  for i=1:size(headername,1)
    text(xc,yc,headername(i,:),'FontWeight','bold','FontSize',fontsize);
    yc=yc-dy*fontsize/12;
  end
end

for i=1:size(varname,1)
	%strip off extra spaces
	j=size(varname,2);
	while (varname(i,j)==' ')
		j=j-1;
	end
	stri=varname(i,1:j);
	
	if (stri(1)=='#')
		text(xc,yc-dy*(i-1),stri(2:length(stri)),'FontName','Symbol');
	else
		text(xc,yc-dy*(i-1),stri);
	end
	text(xc+size(varname,2)*spx,yc-dy*(i-1),sprintf('= %5.2f',var(i)));
end








