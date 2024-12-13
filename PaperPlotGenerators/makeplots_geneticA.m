clear;

f = gcf;
f.Position = [50 50 800 600];

rootdir = '..\QAPdata_combined\';
model = load([rootdir 'model.mat']);
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

outputdir = '.\output_evo\';

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

truenames = ["stf60es2", "sko72", "xran70A1", "term75_4", "stf60er1", "dre72", "stf100ep3", "hyp64_3" ];
bangnames = ["stf60es2", "sko72", "xran70A1", "term75!4", "stf60er1", "dre72", "stf100ep3", "hyp64!3" ];
undernames = ["stf60es2", "sko72", "xran70A1", "term75-4", "stf60er1", "dre72", "stf100ep3", "hyp64-3" ];

sko72 = find(contains(model.data.instlabels,'sko72'));
hyp64 = find(contains(model.data.instlabels,'hyp64_3'));
stf60s = find(contains(model.data.instlabels,'stf60es2'));
stf60r = find(contains(model.data.instlabels,'stf60er1'));
dre72 = find(contains(model.data.instlabels,'dre72'));
term75 = find(contains(model.data.instlabels,'term75_4'));
stf100p = find(contains(model.data.instlabels,'stf100ep3'));
%lipa70b = find(contains(model.data.instlabels,'lipa70b'));
xran70A1 = find(contains(model.data.instlabels,'xran70A1'));

sources = [stf60s;sko72;xran70A1;term75;stf60r;dre72;stf100p;hyp64];

sourcepoints = [model.pilot.Z(sources,:)];

%targets = model.cloist.Zedge(1:end,:);
% Put additional targets inside any large gaps
%targets = fillpath(targets,1.2);

% targets = [-1.5,2; 
%         1,3;
%         2.5,1.5;
%         3, 0.5;
%         3, -1;
%         1.5,-2;
%         0,-2.5;
%         -1.5,-2;
%         -2.75,-1;
%         -2.75,0.5];

targets = [-2, 1.5;
    2, 2.5;
    2.75,0.75;
    2.75, -0.75;
    1.75, -2;
    0.25,-2.5;
    -1.5, -2;
    -2.75,-0.5];

symbols = ['d','s','o','^'];

results = zeros(length(sources),8,2);

newdir = '..\QAP_evoonly\';
fcfeatfile = [newdir 'metadata_evoflow.csv'];
fcfeattable = readtable(fcfeatfile);
fcfeat = table2array(fcfeattable(:,2:end-3));
fcalg = table2array(fcfeattable(:,end-2:end-1));
fcsuppfile = [newdir 'suppdata_evoflow.csv'];
fcsupptable = readtable(fcsuppfile);
fcsubS = categorical(fcsupptable.subsource);

nfcfeat = fcfeat;
if model.opts.auto.preproc && model.opts.bound.flag
    himask = bsxfun(@gt,nfcfeat,model.prelim.hibound);
    lomask = bsxfun(@lt,nfcfeat,model.prelim.lobound);
    nfcfeat = nfcfeat.*~(himask | lomask) + bsxfun(@times,himask,model.prelim.hibound) + ...
                                bsxfun(@times,lomask,model.prelim.lobound);
end

if model.opts.auto.preproc && model.opts.norm.flag
    nfcfeat = bsxfun(@minus,nfcfeat,model.prelim.minX)+1;
    nfcfeat = bsxfun(@rdivide,bsxfun(@power,nfcfeat,model.prelim.lambdaX)-1,model.prelim.lambdaX);
    nfcfeat = bsxfun(@rdivide,bsxfun(@minus,nfcfeat,model.prelim.muX),model.prelim.sigmaX);
end

fcfeatsel = nfcfeat(:,model.featsel.idx);
evolvedZ = fcfeatsel*model.pilot.A';
evolvednames = string(fcfeattable.instances);

for i = 1:length(sources)
    for j = 1:8
        iname = strcat("evoflow_",bangnames(i),"_",num2str(j),"_1");
        iindex = find(iname == evolvednames);
        if isempty(iindex)
            error(iname);
        else
            results(i,j,1) = evolvedZ(iindex,1);
            results(i,j,2) = evolvedZ(iindex,2);
        end
    end
end

% for i = 1:length(sources)
%     for j = 1:8
%         results(i,j,1) = sourcepoints(i,1) + 2*rand - 1;
%         results(i,j,2) = sourcepoints(i,2) + 2*rand - 1;
%     end
% end

i = 1;
nperplot = 4;
nplots = ceil(length(sourcepoints) / nperplot);

xr = results(:,:,1);
xr = xr(:);
yr = results(:,:,2);
yr = yr(:);

minx = min([model.pilot.Z(:,1); targets(:,1); xr ]);
maxx = max([model.pilot.Z(:,1); targets(:,1); xr ]);
miny = min([model.pilot.Z(:,2); targets(:,2); yr ]);
maxy = max([model.pilot.Z(:,2); targets(:,2); yr ]);

for p = 1:nplots

    figure(p)
    f = gcf;
    f.Position = [50 50 800 600];
        
    scatter(model.pilot.Z(:,1), model.pilot.Z(:,2), 10, [0.94 0.94 0.94], 'filled');
    hold on

    for k = 1:size(targets,1)
        thand = scatter(targets(k,1), targets(k,2), 50, [0 0 0], 'p', 'filled');
    end
    
    start = nperplot*(p-1);
    nhere = min(nperplot, length(sourcepoints)-start);

    cmap = @copper;
    clrs = flipud(cmap(nhere+1));
    clrs = clrs(2:end,:);
    clrs = [clrs(1:2:nhere,:); clrs(2:2:nhere,:)];

    hand = [];

    for j = 1:nhere
        i = start + j;
        for k = 1:8
            plot([sourcepoints(i,1),results(i,k,1)],[sourcepoints(i,2),results(i,k,2)], 'LineWidth', 1', 'Color', clrs(j,:));
        end
    end

    for j = 1:nhere
        i = start + j;
        for k = 1:8
            hand(j) = scatter(results(i,k,1), results(i,k,2), 50, symbols(j), 'MarkerEdgeColor', clrs(j,:), 'MarkerFaceColor', [1,1,1], 'LineWidth',1.5);
        end
    end

    for j = 1:nhere
        i = start + j;
        scatter(sourcepoints(i,1), sourcepoints(i,2), 150, symbols(j), 'MarkerEdgeColor', clrs(j,:), 'MarkerFaceColor', [1,1,1], 'LineWidth',2);
    end

    xlim([minx - 0.2, maxx + 0.2]);
    ylim([miny - 0.2, maxy + 0.2]);

    legend([thand,hand],["Targets",undernames(start+1:start+nhere)],"Location","southoutside",'NumColumns', nhere+1);
        
    hold off

    title("Evolved instances")

    print(gcf,'-dpng',strcat(outputdir,'evoonly_spider_',num2str(p),'.png'));
    print(gcf,'-depsc',strcat(outputdir,'evoonly_spider_',num2str(p),'.eps'));

end

