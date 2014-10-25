%% function [im_rgb] = getImage(datadir,fname,opts)
% 
% INPUTS
%  datadir       - directory of H&E images;
%                  later on will be WSI and then will be broken down 
%  fname         - file name
%  opts          - parameter settings (see setEnvironment_input)
%               file type of svs or tile
%
% OUTPUTS
%  im_rgb        - NxMx3 query image
% 

% -------------------------------------------------------------------------
% HE tumor object segmentation 
% Luong Nguyen, lun5@pitt.edu
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------
function [im_rgb] = getImage(datadir,fname,opts)
if strcmp(opts.fileType,'whole slide')
    disp('currently cannnot read this file with JPEG2000 compression');
    im_rgb = [];
elseif opts.fileType = 'tiled'
    im_rgb = imread(fullfile(datadir, fname));
else
    error('Wrong input file type. Need to be either whole slide or tiled');    
end

%% Need to break each image into smaller squares, how did Virginia do this
% right now I am just looking at small images of size 100-500


end