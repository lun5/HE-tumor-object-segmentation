datadir = 'T:\HE_Tissue-Image(Luong)\TissueImages';
if ~ exist(datadir,'dir')
    datadir = '/Users/lun5/Research/color_deconvolution/TissueImages/';
end

resultdir = fullfile(pwd,'test_images');
svs_fname = 'tp10-867-1';
wsi_get_objects(svs_fname, datadir, resultdir); 