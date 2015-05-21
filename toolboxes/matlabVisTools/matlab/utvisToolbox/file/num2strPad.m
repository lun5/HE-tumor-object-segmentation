function s=num2strPad(n, pad)
%
% s=num2strPad(n, pad)
%

  if ~exist('pad')
     pad =1;
  end	
  s = num2str(n);
  while (length(s) < pad)
    s = ['0' s];
  end

