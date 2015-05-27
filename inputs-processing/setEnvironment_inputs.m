%% function [opts_input] = setEnvironment_inputs
% set parameter options for getting input images
% 
% INPUTS
%  fileType - svs or already tiled image 
%  tilesize - 
%  type - specifies which parameter set to use 
%          e.g., can take on values 'speedy' or 'accurate' 
%          feel free to define your custom types at end of this function
%  parallelize option
%  
% OUTPUTS
%  opts_input - selected parameters
% 
% -------------------------------------------------------------------------
% HE segmentation toolbox
% Luong Nguyen, 2014 [lun5@pitt.edu]
% Please email me if you find bugs, or have suggestions or questions
% -------------------------------------------------------------------------

function [opts_input] = setEnvironment_inputs
    %% file type
    opts_input.fileType = 'tiled' ; % default value
    
    %% size of images to process
    opts_input.tileSize = 300;
    
   
end