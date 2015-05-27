function [p] = learnP_A_Blocal(F, opts)
   numSamples = size(F,1);
   F_val = F(randperm(numSamples,500),:);
   p = kde(F','lcv',[],'e');
   f = @(bw,p,F_val) nLOO_LL_anis(bw,p,F_val);
   fminsearch_opts.Display = 'off';%'iter';
   fminsearch_opts.MaxIter = 20;
   reg_min = opts.kde.min_bw; % this regularizes for things like perceptual discriminability, that do not show up in the likelihood fit
                              %  reduces the impact of
                              %  outlier channels and noise
   reg_max = opts.kde.max_bw;
           
   for i=1:2 % for some reason repeatedly running fminsearch continues to improve the objective
        bw = getBW(p,1);
        bw_star = fminsearch(@(bw) f(bw,p,F_val), bw(1:size(bw,1)/2), fminsearch_opts);
        bw_star = cat(1,bw_star,bw_star);
        adjustBW(p,min(max(bw_star,reg_min),reg_max));
   end
end

function [H] = nLOO_LL_anis(bw,p,F_val)
    bw = cat(1,bw,bw);
    adjustBW(p,bw);
    H = -evalAvgLogL(p,F_val');
end