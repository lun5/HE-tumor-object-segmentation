%% matching brightness by matching moments
function source_eq = matchingMoments(source_im, target_im,which_feature) % to be added saturated indices
    %% normalize the range of brightness to between 0 and 1
    target_im = target_im(:);
    source_im = source_im(:);
    if strcmp(which_feature,'brightness opp')|| strcmp(which_feature,'saturation opp')
        % statistics of target
        target_mean = mean(target_im);
        target_std = std(target_im);
        target_skewness = skewness(target_im);
        target_kurtosis = kurtosis(target_im);

        % statistics of source
        source_mean = mean(source_im);
        source_std = std(source_im);
        source_skewness = skewness(source_im);
        source_kurtosis = kurtosis(source_im);

        % fix the mean and variance
        source_eq = (source_im - source_mean)*target_std/source_std + target_mean;
        % fix the range
        source_eq(source_eq <= min(target_im)) = min(target_im);
        source_eq(source_eq >= max(target_im)) = max(target_im);
        % fix skewness 
        [source_eq, ~] = modskew(source_eq,target_skewness);
        % fix kurtosis
        [source_eq, ~] = modkurt(source_eq,target_kurtosis);
        source_eq_skewness = skewness(source_eq);
        source_eq_kurtosis = kurtosis(source_eq);
        %% statistics after equalization
        sprintf('Feature %s \n',which_feature)
        sprintf('Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nMatched Source %.4f %.4f %.4f %.4f %.4f %.4f\n',...
        source_mean, source_std, min(source_im),...
        max(source_im), source_skewness, source_kurtosis,...
        target_mean, target_std, min(target_im),...
        max(target_im), target_skewness, target_kurtosis,...
        mean(source_eq), std(source_eq),...
        min(source_eq), max(source_eq),...
        source_eq_skewness, source_eq_kurtosis)

    elseif strcmp(which_feature,'hue opp')
        % mixture of univariate von Mises distributions
        X_cart_target = [cos(target_im) sin(target_im)];
        %% Call the function
        numClusters = 3; opts_mixture.noise = 1;
        [ mu_hat_polar_target,~, ~,posterior_probs_target, ~] =...
           moVM(X_cart_target,numClusters, opts_mixture); 
        [~, indx_membership_target] = max(posterior_probs_target,[],2); % 4 is the uniform noise
        
        X_cart_source = [cos(source_im) sin(source_im)];
        %% Call the function
        opts_mixture.noise = 1;
        [ mu_hat_polar_source,~, ~,posterior_probs_source, ~] =...
           moVM(X_cart_source,numClusters, opts_mixture); % no noise version
        [~, indx_membership_source] = max(posterior_probs_source,[],2); % 4 is the uniform noise
        source_eq = zeros(size(source_im));
        %% normalize for each cluster
        for cl = 1:(numClusters)% + 1)
            target_im_cl = target_im(indx_membership_target == cl)';
            source_im_cl = source_im(indx_membership_source == cl)';
            % statistics of target
            target_mean = mu_hat_polar_target(cl); %circ_mean(target_im_cl,[],2);
            [~,target_std] = circ_std(target_im_cl,[],[],2);
            
            % statistics of source
            source_mean = mu_hat_polar_source(cl); %circ_mean(source_im_cl,[],2);
            [~,source_std] = circ_std(source_im_cl,[],[],2);
            %% center to 0
            source_eq_cl = source_im - source_mean; %all of the image, not just the cluster of interest
            source_eq_cl(source_eq_cl < - pi) = source_eq_cl(source_eq_cl < - pi)+ 2*pi;
            source_eq_cl(source_eq_cl > pi) = source_eq_cl(source_eq_cl > pi)- 2*pi;

            %% spread source to have same std as target
            source_eq_cl = source_eq_cl * target_std/source_std;

            %% fix mean to be equal target mean
            source_eq_cl = source_eq_cl + target_mean;
            source_eq_cl(source_eq_cl < - pi) = source_eq_cl(source_eq_cl < - pi)+ 2*pi;
            source_eq_cl(source_eq_cl > pi) = source_eq_cl(source_eq_cl > pi)- 2*pi;
            source_eq_cl = source_eq_cl.*posterior_probs_source(:,cl);
            %% put it back together
            source_eq = source_eq + source_eq_cl;
            
            %sprintf('Feature %s Cluster number %d\n',which_feature, cl)
            %sprintf('Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nMatched Source %.4f %.4f %.4f %.4f %.4f %.4f\n',...
            %source_mean, source_std, min(source_im),...
            %max(source_im), source_skewness, source_kurtosis,...
            %target_mean, target_std, min(target_im),...
            %max(target_im), target_skewness, target_kurtosis,...
            %circ_mean(source_eq_cl'), circ_std(source_eq_cl'),...
            %min(source_eq_cl), max(source_eq_cl),...
            %circ_skewness(source_eq_cl'),circ_kurtosis(source_eq_cl'))
        end
        
        if opts_mixture.noise
            source_eq = source_eq + source_im.*posterior_probs_source(:,numClusters + 1);           
        end
        
        source_eq(source_eq < - pi) = source_eq(source_eq < - pi)+ 2*pi;
        source_eq(source_eq > pi) = source_eq(source_eq > pi)- 2*pi;

        
        
        
    end
    
    %% histogram
    A = {source_im,target_im,source_eq};
    legendstr = {'source','target','matched source'};
    bincol = [0.8,0.8,0.8];
    figure; nhist(A,'legend',legendstr,'color',bincol,'xlabel',[which_feature ' values'],...
    'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');
end