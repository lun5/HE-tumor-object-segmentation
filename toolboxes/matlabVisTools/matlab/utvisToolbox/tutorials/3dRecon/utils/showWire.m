function lines = showWire(v, f, h)
  if nargin < 3
    h =1;
  end
  
  figure(h); clf;
  plot3(v(:,1), v(:,2), v(:,3), '.b');
  hold on;
  
  lines = {};
  for k = 1:length(f)
    id = f{k};
    id = [id id(1)];
    lines{k} = v(id ,:);
    l = lines{k};
    plot3(l(:,1), l(:,2), l(:,3), 'g');
  end
  axis vis3d; axis square; axis equal;
