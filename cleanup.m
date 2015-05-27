% this is a script to visualize the code. The purpose is to clean up the
% code but it might also make a mess out of it. 
GitHubdir = 'C:\Users\luong_nguyen\Documents\GitHub\';
addpath(genpath(fullfile(GitHubdir,'plot_subfun20150311')));
addpath(genpath(fullfile(GitHubdir,'HE-tumor-object-segmentation')));

data = plot_subfun('masterscript','-extshow','-hide','parsave');
for iter_fun = 1:length(data.external)
    figure;
    plot_subfun(data.external{iter_fun},'-extshow','-hide','parsave');
end
plot_subfun('masterscript','-unusedonly','-extlist');
