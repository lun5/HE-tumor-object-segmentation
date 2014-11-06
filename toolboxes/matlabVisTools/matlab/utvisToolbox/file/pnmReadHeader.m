function [pars, type]= pnmReadHeader(fid)
% [pars, type]= pnmReadHeader(fid)
% Read and parse a pgm style image file header
% For a successful parse, return:
%  pars = [xdim ydim maxval]
%       xdim, ydim = number of columns, number of rows
%       maxval = maximum pixel value
% and type:
%  "P1" = ascii bitmap, "P2" = ascii greymap,
%  "P3" = ascii pixmap, "P4" = raw bitmap, 
%  "P5" = raw greymap, "P6" = raw colour pixmap
%  "PL", "FP" = float file, ieee-le format
%  "PB", "FU" = float file, ieee-be format
% For an unsuccessful parse, returns 
%  pars=-1
% Written by A. Jepson, 9/99

pars = -1;
type = 'Unknown';

%%% First line contains type string:
TheLine = fgetl(fid);
szLine = size(TheLine);
endLine = szLine(2);

if (endLine < 2)
  fprintf(1, ['Unrecognized PNM file type\n']);
  pars = -1;
  return;
end

type = TheLine(1:2);
ok = 0;
if (type(1) == 'P')
   if (type(2) == '1' | type(2) == '2' | ...
       type(2) == '3' | type(2) == '4' | ...
       type(2) == '5' | type(2) == '6' | ...
       type(2) == 'B' | type(2) == 'L')
     ok = 1;
   else
     fprintf(1, ['Unrecognized PNM file type: ' type '\n']);
   end
end
if (type(1) == 'F')      
   if (type(2) == 'P' | type(2) == 'U')
     ok = 1;
   else
     fprintf(1, ['Unrecognized PNM file type: ' type '\n']);
   end
end
  
if ~ok
   pars = -1;
   return;
end

%% First two characters of first line have been understood.
current = 3; % index of current character in TheLine

parIndex=1;  % index for next parameter to be read
while(parIndex < 4) 
  while (current > endLine)
    TheLine = fgetl(fid);
    if (TheLine == -1)
      fprintf(1, 'Unexpected EOF\n');
      pars = -1;
      return;
    end
    szLine = size(TheLine);
    endLine = szLine(2);
    current = 1;
  end
  
  [token, count, errmsg, nextindex] = ...
       sscanf(TheLine(current:endLine),'%s',1);
  nextindex = nextindex+current-1;
  
  if (count==0)
    if (nextindex > endLine) % Read to end of line
     current = nextindex;
    else % Problem
     pars = -1;
     fprintf(1, 'Unexpected EOF\n');
     return;
    end
  else  % Found a token
  
   if token(1) == '#'
    current = endLine+1;
   else
    [pars(parIndex), count, errmsg, nextindex] = ...
       sscanf(TheLine(current: endLine), '%d', 1);
    if ~(count==1)
      fprintf(1,'Confused reading pgm header\n');
      pars=-1;
      return;
    end
    parIndex = parIndex+1;
    current = current+nextindex-1;
   end

  end
end
