function rgb=name2rgb(name)
%rgb=name2rgb(name)
%Returns rgb values for color names:
%  'y' yellow   [1   1   0]
%  'm' magenta  [1   0   1]
%  'c' cyan     [0   1   1]
%  'r' red      [1   0   0]
%  'g' green    [0   1   0]
%  'b' blue     [0   0   1]
%  'w' white    [1   1   1]
%  'k' black    [0   0   0]

%4/12/96  gmb   wrote it.
collist='ymcrgbwk';

ColorOrder=[
     1     1     0
     1     0     1
     0     1     1
     1     0     0
     0     1     0
     0     0     1
     1     1     1
     0     0     0];

rgb=ColorOrder(find(collist==name),:);

