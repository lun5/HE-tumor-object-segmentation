function ChangeColors(fromcol,tocol,handles)
%ChangeColors(fromcol,tocol,[handles])
%  Changes colors on a graph from 'fromcol' to 'tocol'.
%  'fromcol' and 'tocol' can be either rgb arrays, such as [1,0,0' 
%   or color name strings such as 'y' or 'b'.
%  ChangeColors works recursively through the current graph (gca)
%  and changes the colors of all children with types 'patches', 'lines', 
%  and 'text'.
%
%Example: 
%  ChangeColors('y','b')
%  Changes all yellow lines, patches and text on the graph to blue
%
%  ChangeColors('y',[0.5,0.5,0.5])
%  Changes all yellow lines, patches and text on the graph to grey

%REQUIRES:
%  name2rgb.m, strcomp.m

%4/121/96 gmb  wrote it.

if (nargin<3)
	handles=gcf;
end

if length(fromcol)==1
	fromcol=name2rgb(fromcol);
end

if length(tocol)==1
	tocol=name2rgb(tocol);
end

for i=1:length(handles)
  %some tricky recursive programming here...
  children=get(handles(i),'Children');
  if children~=[];
	ChangeColors(fromcol,tocol,children);
  end

  type=get(handles(i),'Type');
  %disp(type);
  if strcomp(type,'line')==1
	  col=get(handles(i),'Color');
	  if strcomp(col,fromcol)==1	
		col=tocol;
		set(handles(i),'Color',col);
	  end
   end
  if strcomp(type,'text')==1
	  col=get(handles(i),'Color');
	  if strcomp(col,fromcol)==1	
		col=tocol;
		set(handles(i),'Color',col);
	  end
   end
  if strcomp(type,'patch')==1
	  col=get(handles(i),'FaceColor');
	  if strcomp(col,fromcol)==1	
		col=tocol;
		set(handles(i),'FaceColor',col);
	  end
   end
end




