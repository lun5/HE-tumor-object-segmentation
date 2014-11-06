function setFontSize(size,handles)
%
% setFontSize(size,[handles])
%
% djh, 1/24/98

if ~exist('handles')
  handles=gcf;
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
end


