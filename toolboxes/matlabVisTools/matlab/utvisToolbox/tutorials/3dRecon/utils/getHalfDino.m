function [v, f] = getHalfDino(CHATTY)
%% [verts, faces] = getHalfDino(CHATTY)
%% If CHATTY (default FALSE) then do display
FALSE = 0==1;
TRUE = ~FALSE;

if nargin == 0
  CHATTY = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read and display Dino data set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dinoData;

if CHATTY
  lines = showWire(verts,faces,1);
  title('Full Dino Data Set');
  fprintf(2,'Rotate this figure.\n');
  fprintf(2,'Press any key to continue...');
  pause; fprintf(2,'ok\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract one side of Dino in principal coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Dino's principal axes
mn = sum(verts,1)/size(verts,1);
[U S V] = svd(verts-repmat(mn, size(verts,1), 1), 0); S = diag(S);

%% Expect third axis to be the lateral direction.
%% Remove vertices with z < 0
id = U(:,2) < -0.001; 

%% Compute vertices in principal axes. 
v = U(~id, :)*diag(S);

%% Relabel the remaining vertices
relabel = zeros(size(verts,1),1);
relabel(~id) = 1:sum(~id);

%% Write the faces in terms of the relabelled vertices
f={}; j=0;
for k = 1:length(faces)
  if ~any(id(faces{k}))
    j = j+1;
    f{j} = relabel(faces{k})';
  end
end

if CHATTY
  l = showWire(v,f,1);
  title('One Sided Dino Data Set');
  fprintf(2,'Rotate this figure.\n');
  fprintf(2,'Press any key to continue...');
  pause; fprintf(2,'ok\n');

  %%%%% Show Dino as a surface plot
  figure(3); clf;
  for k = 1:length(f)
    vf = v(f{k},:);
    patch(vf(:,1), vf(:,2), vf(:,3), vf(:,2));
  end
  set(gca,'YDir', 'reverse');
  axis vis3d; axis square; axis equal;
  title('One Sided Dino Data Set');
  fprintf(2,'Rotate this figure.\n');
  fprintf(2,'Press any key to continue...');
  pause;
  fprintf(2,'ok\n');

end
return;