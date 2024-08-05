% Project feature data for a new set of instances (with algorithm data)
% onto the existing 2-D instance space.

rootdir = '.\QAPdata\';
model = load('.\QAPdata\model.mat');

suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

newmetafile = [rootdir 'metadata_recomb.csv'];
newmetatable = readtable(newmetafile);

varlabels = newmetatable.Properties.VariableNames;
isname = strcmpi(varlabels,'instances');
isfeat = strncmpi(varlabels,'feature_',8);
isalgo = strncmpi(varlabels,'algo_',5);
issource = strcmpi(varlabels,'source');
newinstlabels = newmetatable{:,isname};
if isnumeric(model.data.instlabels)
    newinstlabels = num2cell(model.data.instlabels);
    newinstlabels = cellfun(@(x) num2str(x),model.data.instlabels,'UniformOutput',false);
end
if any(issource)
    newsources = categorical(newmetatable{:,issource});
end
newfeat = newmetatable{:,isfeat};
newalg = newmetatable{:,isalgo};


%newfeat = table2array(newfeattable(:,2:end));
%newfeatraw = newfeat;

if model.opts.auto.preproc && model.opts.bound.flag
    himask = bsxfun(@gt,newfeat,model.prelim.hibound);
    lomask = bsxfun(@lt,newfeat,model.prelim.lobound);
    newfeat = newfeat.*~(himask | lomask) + bsxfun(@times,himask,model.prelim.hibound) + ...
                                bsxfun(@times,lomask,model.prelim.lobound);
end

if model.opts.auto.preproc && model.opts.norm.flag
    newfeat = bsxfun(@minus,newfeat,model.prelim.minX)+1;
    newfeat = bsxfun(@rdivide,bsxfun(@power,newfeat,model.prelim.lambdaX)-1,model.prelim.lambdaX);
    newfeat = bsxfun(@rdivide,bsxfun(@minus,newfeat,model.prelim.muX),model.prelim.sigmaX);
end

newfeatsel = newfeat(:,model.featsel.idx);
newfeatZ = newfeatsel*model.pilot.A';

%figure
clf
%set(gcf, 'Position',  [0, 100, 800, 800])
Z = model.pilot.Z;
scatter(Z(:,1), Z(:,2),10,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
hold on
scatter(newfeatZ(:,1), newfeatZ(:,2),120,'Marker','x','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
hold off
