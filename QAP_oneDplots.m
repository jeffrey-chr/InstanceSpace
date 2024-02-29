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

for alg = 1:length(model.data.algolabels)
    for i = 2:size(unft,2)
        tmp = unft.Properties.VariableNames(i);
        scatter(unft.(tmp{1}), model.data.Yraw(mnI,alg),'blue','.');
        title([tmp{1}(9:end) ' vs ' model.data.algolabels{alg}] );
        %xlim([0, max(unft.feature_Size) + 10])
        xlabel('Feature');
        ylabel(model.data.algolabels{alg});
        %print(gcf,'-depsc',['./oneDplots/size/size_' tmp{1}(9:end)]);
        print(gcf,'-dpng',['./oneDplots/' model.data.algolabels{alg} '/' model.data.algolabels{alg} '_' tmp{1}(9:end)]);
    end
end