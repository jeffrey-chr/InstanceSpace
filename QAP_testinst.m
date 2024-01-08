% Project feature data for a new set of instances (without algorithm data)
% onto the existing 2-D instance space.

rootdir = '.\QAPdata\';
model = load('.\QAPdata\model.mat');

suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

newfeatfile = [rootdir 'featdata.csv'];
newfeattable = readtable(newfeatfile);
newfeat = table2array(newfeattable(:,2:end));
newfeatraw = newfeat;

if model.opts.auto.preproc && model.opts.bound.flag
    himask = bsxfun(@gt,newfeat,model.bound.hibound);
    lomask = bsxfun(@lt,newfeat,model.bound.lobound);
    newfeat = newfeat.*~(himask | lomask) + bsxfun(@times,himask,model.bound.hibound) + ...
                                bsxfun(@times,lomask,model.bound.lobound);
end

if model.opts.auto.preproc && model.opts.norm.flag
    newfeat = bsxfun(@minus,newfeat,model.norm.minX)+1;
    newfeat = bsxfun(@rdivide,bsxfun(@power,newfeat,model.norm.lambdaX)-1,model.norm.lambdaX);
    newfeat = bsxfun(@rdivide,bsxfun(@minus,newfeat,model.norm.muX),model.norm.sigmaX);
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
