%% function [opts_clustering] = setEnvironment_clustering
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

function [opts_clustering] = setEnvironment_clustering
   %% globalization (from affinity to boundaries)              used in getE:
    opts_clustering.globalization_method = 'spectral_clustering';          % how to go from affinty to boundaries? (spectral clustering is only method currently supported)
    opts_clustering.spectral_clustering.approximate = false;                % use the DNcuts approximation from Arbelaez et al. CVPR 2014? (was not included in our published paper)
    opts_clustering.spectral_clustering.nvec = 100;                        % how many eigenvectors to use
    
    opts_clustering.display_progress = true;
    % post-processing
    opts_clustering.border_suppress = 1;                        % get rid of boundaries that align with image borders and are right next to the borders?
                                                                %  (this helps on images that have false boundaries near borders (like some in BSDS); this kind of suppression
                                                                %   is common in other boundary detection algorithms such as Structured Edges (Dollar & Zitnick 2013) and
                                                                %   Sketch Tokens (Lim et al. 2013))
                                     

end