clear;

f = gcf;
f.Position = [50 750 800 600];

cmap = @copper;

alteredscriptfcn;

rootdir = '..\QAPdata\';
model = load([rootdir 'model.mat']);
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

newdir = '..\QAP_separateexp\';
fcfeatfile = [newdir 'metadata_fc.csv'];
fcfeattable = readtable(fcfeatfile);
fcfeat = table2array(fcfeattable(:,2:end-3));
fcalg = table2array(fcfeattable(:,end-2:end-1));
fcsuppfile = [newdir 'suppdata_fc.csv'];
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
fcfeatZ = fcfeatsel*model.pilot.A';


recombfeatfile = [newdir 'metadata_recomb.csv'];
recombfeattable = readtable(recombfeatfile);
recombfeat = table2array(recombfeattable(:,2:end-3));
recombalg = table2array(recombfeattable(:,end-2:end-1));
recombsuppfile = [newdir 'suppdata_recomb.csv'];
recombsupptable = readtable(recombsuppfile);
rcsubS = categorical(recombsupptable.subsource);

nrcfeat = recombfeat;
if model.opts.auto.preproc && model.opts.bound.flag
    himask = bsxfun(@gt,nrcfeat,model.prelim.hibound);
    lomask = bsxfun(@lt,nrcfeat,model.prelim.lobound);
    nrcfeat = nrcfeat.*~(himask | lomask) + bsxfun(@times,himask,model.prelim.hibound) + ...
                                bsxfun(@times,lomask,model.prelim.lobound);
end

if model.opts.auto.preproc && model.opts.norm.flag
    nrcfeat = bsxfun(@minus,nrcfeat,model.prelim.minX)+1;
    nrcfeat = bsxfun(@rdivide,bsxfun(@power,nrcfeat,model.prelim.lambdaX)-1,model.prelim.lambdaX);
    nrcfeat = bsxfun(@rdivide,bsxfun(@minus,nrcfeat,model.prelim.muX),model.prelim.sigmaX);
end

rcfeatsel = nrcfeat(:,model.featsel.idx);
rcfeatZ = rcfeatsel*model.pilot.A';

outputdir = '.\output\';

% flowcluster plot

pfcZ = [model.pilot.Z; fcfeatZ];
subs = [supp.subsource; fcsupptable.subsource];

fcsources = repmat([""], length(subs), 1);
for i = 1:length(fcsources)
    if startsWith(subs{i},"flowcluster-dhyper-fcycle")
        fcsources(i) = "Hypercube x 3-cycle";
    elseif startsWith(subs{i},"flowcluster-dhyper-ftree")
        fcsources(i) = "Hypercube x Tree";
    elseif startsWith(subs{i},"flowcluster-dhyper-fsquare")
        fcsources(i) = "Hypercube x Square";
    elseif startsWith(subs{i},"flowcluster-ddrez-fcycle")
        fcsources(i) = "Drexx x 3-cycle";
    elseif startsWith(subs{i},"flowcluster-ddrez-ftree")
        fcsources(i) = "Drexx x Tree";
    elseif startsWith(subs{i},"flowcluster-ddrez-fsquare")
        fcsources(i) = "Drexx x Square";
    else
        fcsources(i) = "";
    end
end
fcsourcescat = categorical(fcsources);
typs = {"Hypercube x 3-cycle", "Hypercube x Tree", "Hypercube x Square","Drexx x 3-cycle", "Drexx x Tree", "Drexx x Square"};

drawSources(pfcZ, fcsourcescat, cmap, typs);
title('New Flow-Cluster Instances')
print(gcf,'-dpng',[outputdir 'iniplus_flowcluster.png']);
print(gcf,'-depsc',[outputdir 'iniplus_flowcluster.eps']);

% recomb plot

% flowcluster plot

pfcZ = [model.pilot.Z; rcfeatZ];
subs = [supp.subsource; recombsupptable.subsource];

rcsources = repmat([""], length(subs), 1);
for i = 1:length(rcsources)
    if startsWith(subs{i},"recomb")
        rcsources(i) = "Recombined";
    else
        rcsources(i) = "";
    end
end
rcsourcescat = categorical(rcsources);
typs = {"Recombined"};

drawSources(pfcZ, rcsourcescat, cmap, typs);
title('New Recombined Instances')
print(gcf,'-dpng',[outputdir 'iniplus_recomb.png']);
print(gcf,'-depsc',[outputdir 'iniplus_recomb.eps']);
