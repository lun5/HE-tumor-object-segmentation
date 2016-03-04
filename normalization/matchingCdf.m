%% matching brightness by cdf
function source_eq_nw = matchingCdf(source_im, target_im,opts) % to be added saturated indices
    %% normalize the range of brightness to between 0 and 1
    target_im = target_im(:);
    source_im = source_im(:);
    indx_white_source = opts.source_stats.posterior_probs(:,3);
    indx_white_target = opts.target_stats.posterior_probs(:,3);
    % not white pixels
    target_im_nw = target_im(~indx_white_target);
    source_im_nw = source_im(~indx_white_source);
    % normalize the range of brightness to between 0 and 1
    target_norm = target_im_nw./range(target_im_nw) + min(target_im_nw);
    source_norm = source_im_nw./range(source_im_nw) + min(source_im_nw);
    % count number of elements in each bin 
    binranges = 0:0.01:1;
    bincounts = histc(target_norm,binranges);
    % equalize source's brightness --> target's brightness in [0 1]
    source_norm_eq = histeq(source_norm,bincounts);
    % multiply by the range of target's brightness to complete normalization
    source_eq_nw = (source_norm_eq - min(target_im_nw))* range(target_img);
    
    %% statistics after equalization
    sprintf('Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nSource_cdf %.4f %.4f %.4f %.4f %.4f %.4f\n',...
    mean(source_im_nw), std(source_im_nw), min(source_im_nw),...
    max(source_im_nw), skewness(source_im_nw), kurtosis(source_im_nw),...
    mean(target_im_nw), std(target_im_nw), min(target_brightness),...
    max(target_im_nw), skewness(target_im_nw), kurtosis(target_im_nw),...
    mean(source_eq_nw), std(source_eq_nw),...
    min(source_eq_nw), max(source_eq_nw),...
    skewness(source_eq_nw), kurtosis(source_eq_nw))

    source_eq = source_im;
    source_eq(~indx_white_source) = source_eq_nw;
    %% histogram
    A = {source_im_nw,target_im_nw,source_eq_nw};
    legendstr = {'source','target','matched source'};bincol = [0.8,0.8,0.8];
    figure; nhist(A,'legend',legendstr,'color',bincol,'xlabel',[which_feature ' values'],...
    'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');
end