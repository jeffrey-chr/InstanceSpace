%% Example Title
% Summary of example objective

rootdir = '.\QAPdata\';
model = load('.\QAPdata\model.mat');
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
unftfile = [rootdir 'UNSCALED_featdata.csv'];
unft = readtable(unftfile);

% match instance labels

[names, mnI] = sort(model.data.instlabels);
[~, ufI] = sort(unft.instances);

if ~all(strcmp(model.data.instlabels(mnI), unft.instances(ufI)))
    error("Instance names do not match");
end

for i = 2:size(unft,2)
    tmp = unft.Properties.VariableNames(i);
    scatter(unft.feature_Size, unft.(tmp{1}),'blue','filled');
    title(tmp{1}(9:end));
    xlim([0, max(unft.feature_Size) + 10])
    xlabel('Size');
    ylabel('Feature');
    %print(gcf,'-depsc',['./oneDplots/size/size_' tmp{1}(9:end)]);
    print(gcf,'-dpng',['./oneDplots/size/size_' tmp{1}(9:end)]);
end

for i = 2:size(unft,2)
    tmp = unft.Properties.VariableNames(i);
    scatter(unft.(tmp{1}), model.data.Yraw(mnI,1),'blue','.');
    hold on
    scatter(unft.(tmp{1}), model.data.Yraw(mnI,2),'red','.');
    scatter(unft.(tmp{1}), model.data.Yraw(mnI,3),'green','.');
    hold off
    legend(model.data.algolabels{1},model.data.algolabels{2},model.data.algolabels{3},'Location','eastoutside');
    title(tmp{1}(9:end));
    %xlim([0, max(unft.feature_Size) + 10])
    xlabel('Feature');
    ylabel('Performance');
    %print(gcf,'-depsc',['./oneDplots/size/size_' tmp{1}(9:end)]);
    print(gcf,'-dpng',['./oneDplots/perf/perf_' tmp{1}(9:end)]);
end