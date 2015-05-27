function [mixture_params] = learnMixtureModel(F,opts)
% calculate pmi internal for hue opp
    if opts.model_half_space_only
        [ params,~, prior_probs,init_params] = mixture_of_bivariate_VM(F, 6);
    else
        [ params,~, prior_probs,init_params] = mixture_of_bivariate_VM(F, 9);
    end                       
    mixture_params.params = params;
    mixture_params.prior_probs = prior_probs;
    mixture_params.init_params = init_params;
end