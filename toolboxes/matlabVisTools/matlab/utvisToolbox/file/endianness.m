function resp = endianess();
% return endianess of machine
%

%% Check output format for stdout
[n p format] = fopen(1);
resp = 'UNKNOWN';
if (strcmp(format, 'ieee-le'))
  resp = 'LITTLE';
elseif (strcmp(format, 'ieee-be'))
  resp = 'BIG';
elseif (strcmp(format, 'native')) 
  %% Check for a PCWIN computer type, or an Linux x86 indicator
  sys = computer;
  if (findstr('PCWIN', sys)>0) | (findstr('LNX86', sys) > 0)
    resp = 'LITTLE';
  else
    resp = 'BIG';
  end
end
