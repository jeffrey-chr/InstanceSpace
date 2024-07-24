% Project feature data for a new set of instances (without algorithm data)
% onto the existing 2-D instance space.

rootdir = '.\QAPdata\';
%newdir = './Recomb_data/';
newdir = './Recomb_data/';
model = load('.\QAPdata\model.mat');

suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

newfeatfile = [newdir 'metadata_recomb.csv'];
newfeattable = readtable(newfeatfile);
newfeat = table2array(newfeattable(:,2:end-3));
newalg = table2array(newfeattable(:,end-2:end-1));
newfeatraw = newfeat;

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

BMAbad = (newalg(:,1) > 1);
MMASbad = (newalg(:,2) > 1);
AllGood = (BMAbad + MMASbad == 0);

%figure
clf
%set(gcf, 'Position',  [0, 100, 800, 800])
Z = model.pilot.Z;
scatter(Z(:,1), Z(:,2),10,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
hold on
scatter(newfeatZ(MMASbad,1), newfeatZ(MMASbad,2),120,'Marker','x','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
scatter(newfeatZ(BMAbad,1), newfeatZ(BMAbad,2),120,'Marker','+','MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1]);
scatter(newfeatZ(AllGood,1), newfeatZ(AllGood,2),40,'Marker','*','MarkerEdgeColor',[0.6 0.3 0.6],'MarkerFaceColor',[0.6 0.3 0.6]);

hold off
