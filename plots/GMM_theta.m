measures = theta; %brightness;
measures = double(measures(:));
AIC = zeros(1,3);
BIC = zeros(1,3);
obj = cell(1,3);
options = statset('Display','final','MaxIter',200);
for k = 1:3
    obj{k} = gmdistribution.fit(measures,k,'Options',options);
    AIC(k)= obj{k}.AIC;
    BIC(k) = obj{k}.BIC;
end

[minAIC,numComponents_AIC] = min(AIC);
[minBIC,numComponents_BIC] = min(BIC);
numComponents = max(numComponents_AIC,numComponents_BIC);
numComponents

model = obj{numComponents};
mus = model.mu;
sigmas = model.Sigma;
disp('Means are'); mus(:)
disp('Standard devs are'); sigmas(:)

dataRange = max(measures) - min(measures);
binranges = min(measures):dataRange/50:max(measures);
col_vec = {'r','g','b'};

figure; 
    h = histogram(measures,'DisplayStyle','bar' );
    h.FaceColor = [0.8 .8 .8]; h.BinWidth= dataRange/50;
    h.Normalization = 'probability'; 
    axis([min(measures) max(measures) 0 .15]);
hold on;
y = cell(1,numComponents);
x = binranges;
for i = 1:numComponents
    y{i} = normpdf(x,mus(i),sigmas(i));
    plot(x,y{i}./length(x),col_vec{i},'LineWidth',3); hold on;
end
hold off;

tol = opts.kde.kdtree_tol;
p = kde(measures,0.05,[],'e');
pd = evaluate(p,x(:)',tol);
figure; 
bar(binranges,bincounts/sum(bincounts),'histc');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8 .8 .8],'EdgeColor','w');
xlabel('measures');
hold on
plot(x(:),pd./10,'LineWidth',3);
hold off
