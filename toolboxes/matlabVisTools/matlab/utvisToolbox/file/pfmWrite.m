function [] = pfmWrite(mtx, fname, comment);
% pfmWrite(mtx, filename, comment)
%
% Write a MatLab matrix to a pfm (graylevel image) file.
% comments (optional) a string to be included in the pnm style header.

if ~exist('comment')
  comment = '';
end

[fid,msg] = fopen( fname, 'w' );

if (fid == -1)
  error(msg);
end

%%% First line contains ID string:
%%% "PB" = big endian float file
%%% "PL" = big endian float file

endian = endianness;
if strcmp(endian,'LITTLE')
  fprintf(fid,'PL\n');
elseif strcmp(endian,'BIG')
  fprintf(fid,'PB\n');
else
  error(strcat('Unknown endianness: ', endianness));
end


fprintf(fid,'# MatLab PFMWRITE file, saved %s\n',date);
%% Output comments
if length(comment)>0
  EOL = 10;  % End of line
  if comment(end) ~= EOL
    comment = [comment EOL];  % Add an end of line character to the comment.
  end
  lineEnds = findstr(EOL, comment);
  for j=1:length(lineEnds)
    if j == 1
      line = ['# ' comment(1:lineEnds(j))];
    else
      line = ['# ' comment((lineEnds(j-1)+1):lineEnds(j))];
    end
    fprintf(fid, line);
  end
end

%%% dimensions
fprintf(fid,'%d %d\n',size(mtx,2),size(mtx,1));

%%% Maximum pixel value (meaningless in pfm images, but conforms to pnm
%%% style headers)
fprintf(fid,'  0\n');


count  = fwrite(fid, mtx','float');

fclose(fid);

if (count ~= size(mtx,1)*size(mtx,2))
  fprintf(1,'Warning: File output terminated early!');
end
	  
