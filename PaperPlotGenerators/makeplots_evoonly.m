clear;

close all

f = gcf;
f.Position = [50 750 800 600];

cmap = @copper;

alteredscriptfcn;

rootdir = '..\QAPdata_combined\';
model = load([rootdir 'model.mat']);
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

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
fcfeatZ = fcfeatsel*model.pilot.A';

outputdir = '.\output_evo\';

% evoflow plot

pfcY = [model.data.Yraw; fcalg];
pfcZ = [model.pilot.Z; fcfeatZ];
subs = [supp.subsource; fcsupptable.subsource];

% "EvolvedFlow-stf60er1", "^evoflow.*stf.*er.*$";
%     "EvolvedFlow-stf60es2", "^evoflow.*stf.*es.*$";
%     "EvolvedFlow-hyp64!3", "^evoflow.*hyp.*$";
%     "EvolvedFlow-xran70A1", "^evoflow.*xran.*$";
%     "EvolvedFlow-dre72", "^evoflow.*dre.*$";
%     "EvolvedFlow-sko72", "^evoflow.*sko.*$";
%     "EvolvedFlow-term75!4", "^evoflow.*term.*$";

fcsources = repmat([""], length(subs), 1);
for i = 1:length(fcsources)
    if startsWith(subs{i},"EvolvedFlow-stf60er1")
        fcsources(i) = "stf60er1";
    elseif startsWith(subs{i},"EvolvedFlow-stf60es2")
        fcsources(i) = "stf60es2";
    elseif startsWith(subs{i},"EvolvedFlow-hyp64!3")
        fcsources(i) = "hyp64-3";
    elseif startsWith(subs{i},"EvolvedFlow-xran70A1")
        fcsources(i) = "xran70A1";
    elseif startsWith(subs{i},"EvolvedFlow-dre72")
        fcsources(i) = "dre72";
    elseif startsWith(subs{i},"EvolvedFlow-sko72")
        fcsources(i) = "sko72";
    elseif startsWith(subs{i},"EvolvedFlow-term75!4")
        fcsources(i) = "term75-4";
    elseif startsWith(subs{i},"EvolvedFlow-stf100ep3")
        fcsources(i) = "stf100ep3";
    else
        fcsources(i) = "";
    end
end
fcsourcescat = categorical(fcsources);
typs = {"stf60er1", "stf60es2", "hyp64-3","xran70A1", "dre72", "sko72", "term75-4", "stf100ep3"};

drawSources(pfcZ, fcsourcescat, cmap, typs);
title('New Evolved Instances')
print(gcf,'-dpng',[outputdir 'plusevo.png']);
print(gcf,'-depsc',[outputdir 'plusevo.eps']);

% algorithm performance plots
eligible = (fcsources ~= "");
diffYraw = pfcY(:,2) - pfcY(:,1);

clf

hold on
scatter(pfcZ(~eligible,1), pfcZ(~eligible,2), 5, [0.87 0.87 0.87], 'o', 'filled');
[h1, h2, h3] = drawScatterYraw(pfcZ(eligible,:), diffYraw(eligible,:), "Actual Algorithm Performance", cmap);
perfleg = legend([h1,h2,h3], ["BMA stronger", "Similar performance", "MMAS stronger"], 'Location', 'southoutside', 'NumColumns', 3);
cpos = perfleg.Position;
perfleg.set("Position", [(1 - cpos(3))/2, 0.01, cpos(3), cpos(4)]);
gca().set("Position", gca().Position + [0.02 0.04 -0.04 -0.04]);
title("Algorithm Performance, evolved instances")
%gcf.Position(4) = 650;
%perfleg.Position(2) = 0;

%
print(gcf,'-dpng',[outputdir 'evoonly_realperf.png']);
print(gcf,'-depsc',[outputdir 'evoonly_realperf.eps']);

% recomb plot

% flowcluster plot