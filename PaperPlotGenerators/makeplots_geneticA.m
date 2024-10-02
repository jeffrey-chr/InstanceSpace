clear;

rootdir = '..\QAPdata_combined\';
model = load([rootdir 'model.mat']);
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

outputdir = '.\output\';

nfeats = length(model.data.featlabels);


% newmetafile = [rootdir 'metadata_gen.csv'];
% newmetatable = readtable(newmetafile);
% 
% varlabels = newmetatable.Properties.VariableNames;
% isname = strcmpi(varlabels,'instances');
% isfeat = strncmpi(varlabels,'feature_',8);
% isalgo = strncmpi(varlabels,'algo_',5);
% issource = strcmpi(varlabels,'source');
% newinstlabels = newmetatable{:,isname};
% if isnumeric(model.data.instlabels)
%     newinstlabels = num2cell(model.data.instlabels);
%     newinstlabels = cellfun(@(x) num2str(x),model.data.instlabels,'UniformOutput',false);
% end
% if any(issource)
%     newsources = categorical(newmetatable{:,issource});
% end
% newfeat = newmetatable{:,isfeat};
% newalg = newmetatable{:,isalgo};



%sko49, manhattan on left side
%dre56, drezner on right
%stf60es2, real-like on right
%stf80er3 or stf60er1, real-like on left
%hyp64_3, hypercube in bottom right
%term45_9, terminal in bottom right

sko72 = find(contains(model.data.instlabels,'sko72'));
hyp64 = find(contains(model.data.instlabels,'hyp64_3'));
stf60s = find(contains(model.data.instlabels,'stf60es2'));
stf60r = find(contains(model.data.instlabels,'stf60er5'));
dre72 = find(contains(model.data.instlabels,'dre72'));
term75 = find(contains(model.data.instlabels,'term75_4'));
%lipa70b = find(contains(model.data.instlabels,'lipa70b'));
xran70A1 = find(contains(model.data.instlabels,'xran70A1'));

sources = [sko72;dre72;hyp64;term75;stf60s;stf60r;xran70A1];

sourcepoints = [model.pilot.Z(sources,:)];

targets = model.cloist.Zedge(1:end,:);
% Put additional targets inside any large gaps
targets = fillpath(targets,1.2);

symbols = ['d','s','o','^'];

results = zeros(length(sources),8,2);

for i = 1:length(sources)
    for j = 1:8
        results(i,j,1) = sourcepoints(i,1) + 2*rand - 1;
        results(i,j,2) = sourcepoints(i,2) + 2*rand - 1;
    end
end

i = 1;
nperplot = 4;
nplots = ceil(length(sourcepoints) / nperplot);

for p = 1:nplots

    figure(p)
        
    scatter(model.pilot.Z(:,1), model.pilot.Z(:,2), 10, [0.7 0.7 0.7], 'filled');
    hold on

    for k = 1:size(targets,1)
        scatter(targets(k,1), targets(k,2), 50, [1 0 0], 'p', 'filled');
    end
    
    start = nperplot*(p-1);
    nhere = min(nperplot, length(sourcepoints)-start);

    for j = 1:nhere
        i = start + j;
        for k = 1:8
            plot([sourcepoints(i,1),results(i,k,1)],[sourcepoints(i,2),results(i,k,2)], 'LineWidth', 1', 'Color', [1,0,0]);
        end
    end

    for j = 1:nhere
        i = start + j;
        for k = 1:8
            scatter(results(i,k,1), results(i,k,2), 50, symbols(j), 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1,1,1], 'LineWidth',1.5);
        end
    end

    for j = 1:nhere
        i = start + j;
        scatter(sourcepoints(i,1), sourcepoints(i,2), 150, symbols(j), 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1,1,1], 'LineWidth',2);
    end
        
    hold off

end

