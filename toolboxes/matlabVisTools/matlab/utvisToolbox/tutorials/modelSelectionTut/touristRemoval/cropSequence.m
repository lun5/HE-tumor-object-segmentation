function data = cropSequence(im, cropBox)
% data = cropSequence(im, cropBox)
% 
% Read a sequence of pgm image files and return the k^th cropped
% image in data(:,:,k).
% Input:
%   Filenames:
%     rootName####tailName
%   where #### is the image sequence number consisting
%   of exactly pad digits.  The sequence number is assumed
%   to be padded with zeros on the left, if necessary.
%   
%   seqRange  the range of sequence numbers
% 
%   cropBox = [xmin ymin xmax ymax]
% 
% Output
%   data(:,:,k) corresponds to pixels  xmin:xmax, ymin:ymax from
%               frame number seqRange(k)


 data = [];
 for imNum = 1:size(im,3);
   if size(im,1) > 0
     newData = im(cropBox(2):cropBox(4), cropBox(1):cropBox(3), imNum);
     if size(data,1) == 0
       data = newData;
     else
       data = cat(3, data, newData);
     end
   end
 end
