%% embedding PMI into 2D pca
% Gather all the parameters: prior + von mises parameters

% query all the different tumor/stroma/inflammation pieces in the database
% need to put a limits on sizes of these images
% remember their labels in one array/cell
% run mixture model (before and after normalization) to get the params
% put these parameters as features for each of the piece (how much 
% run PCA on these images
% get the locations on the new 2D coordinates
% use Greg code to plot them on a plane