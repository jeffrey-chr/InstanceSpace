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

scatter(model.pilot.Z(:,1), model.pilot.Z(:,2), 10, [0.4 0.4 0.4], 'filled');

%sko49, manhattan on left side
%dre56, drezner on right
%stf60es2, real-like on right
%stf80er3 or stf60er1, real-like on left
%hyp64_3, hypercube in bottom right
%term45_9, terminal in bottom right

sko49 = find(contains(model.data.instlabels,'sko49'));
hyp64 = find(contains(model.data.instlabels,'hyp64_3'));
stf60s = find(contains(model.data.instlabels,'stf60es2'));
stf60r = find(contains(model.data.instlabels,'stf60er1'));
dre56 = find(contains(model.data.instlabels,'dre56'));
term45 = find(contains(model.data.instlabels,'term45_9'));

sources = [sko49;dre56;hyp64;term45;stf60s;stf60r];

sourcepoints = [model.pilot.Z(sources,:)]

targets = model.cloist.Zedge(1:end,:);
% Put additional targets inside any large gaps
targets = fillpath(targets,1.2);

hold on
scatter(sourcepoints(:,1), sourcepoints(:,2), 150, [1 0 0], 'p', 'filled');

symbols = ["o","^","square",">","diamond"];

for j = 1:5
    i = foo(j);
    scatter(targets(i,1), targets(i,2), 50, [1 0 0], symbols(j), 'filled');
end

hold off

figure

scatter(model.pilot.Z(:,1), model.pilot.Z(:,2), 10, [0.4 0.4 0.4], 'filled');

tab30 = find(contains(model.data.instlabels,'xtab30N1'));

targets = model.cloist.Zedge(1:end,:);
% Put additional targets inside any large gaps
targets = fillpath(targets,1.2);

foo = [4,6,8,10,15];

hold on
scatter(model.pilot.Z(tab30,1), model.pilot.Z(tab30,2), 150, [1 0 0], 'p', 'filled');

symbols = ["o","^","square",">","diamond"];

for j = 1:5
    i = foo(j);
    scatter(targets(i,1), targets(i,2), 50, [1 0 0], symbols(j), 'filled');
end

hold off
