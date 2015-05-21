function setFontSize(size,handles)
%
% setFontSize(size,[handles])
%
% Recursively sets fontSize of xTickLabel and yTickLabel and any
% other text.  Default (if handles is not provided) is to start
% with the current axes (gca).
%
% djh, 1/24/98

if ~exist('handles')
  handles=gca;
end

for i=1:length(handles)
  children=get(handles(i),'Children');
  if ~isempty(children);
    setFontSize(size,children);
  end

  type=get(handles(i),'Type');
  if strcomp(type,'text')
    set(handles(i),'FontSize',size);
  end
  if strcomp(type,'axes')
    set(handles(i),'FontSize',size);
    set(handles(i),'xTickLabel',get(handles(i),'xTickLabel'));
    set(handles(i),'yTickLabel',get(handles(i),'yTickLabel'));
  end
end


