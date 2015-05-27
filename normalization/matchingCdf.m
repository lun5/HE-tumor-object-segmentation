%% matching brightness by cdf
function source_eq = matchingCdf(source_im, target_im) % to be added saturated indices
    %% normalize the range of brightness to between 0 and 1
    target_im = target_im(:);
    source_im = source_im(:);
    % normalize the range of brightness to between 0 and 1
    target_norm = target_im./range(target_im) + min(target_im);
    source_norm = source_im./range(source_im) + min(source_im);
    % count number of elements in each bin 
    binranges = 0:0.01:1;
    bincounts = histc(target_norm,binranges);
    % equalize source's brightness --> target's brightness in [0 1]
    source_norm_eq = histeq(source_norm,bincounts);
    % multiply by the range of target's brightness to complete normalization
    source_eq = (source_norm_eq - min(target_im))* range(target_img);
    
    %% statistics after equalization
    sprintf('Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nSource_cdf %.4f %.4f %.4f %.4f %.4f %.4f\n',...
    mean(source_im), std(source_im), min(source_im),...
    max(source_im), skewness(source_im), kurtosis(source_im),...
    mean(target_im), std(target_im), min(target_brightness),...
    max(target_im), skewness(target_im), kurtosis(target_im),...
    mean(source_eq), std(source_eq),...
    min(source_eq), max(source_eq),...
    skewness(source_eq), kurtosis(source_eq))
    %% histogram
    A = {source_im,target_im,source_eq};
    legendstr = {'source','target','matched source'};bincol = [0.8,0.8,0.8];
    figure; nhist(A,'legend',legendstr,'color',bincol,'xlabel',[which_feature ' values'],...
    'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');
end