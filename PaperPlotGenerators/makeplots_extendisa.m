% Create plots of intitial space 
clear;

rootdir = '..\QAPdata_combined\';
model = load([rootdir 'model.mat']);
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

outputdir = '.\output_extisa\';

nfeats = length(model.data.featlabels);

f = gcf;
f.Position = [400 400 800 600];

cmap = @copper;

alteredscriptfcn

subs = supp.subsource;

% source plots

% qaplib plot
qsources = repmat([""], length(supp.subsource), 1);
for i = 1:length(qsources)
    if contains(supp.subsource{i},"qaplib")
        qsources(i) = "QAPLIB instances";
    else
        qsources(i) = "";
    end
end
qsourcescat = categorical(qsources);
typs = {"QAPLIB instances"};

drawSources(model.pilot.Z, qsourcescat, cmap, typs);
title('QAPLIB Instances')
print(gcf,'-dpng',[outputdir 'extisa_qaplib.png']);
print(gcf,'-depsc',[outputdir 'extisa_qaplib.eps']);


% category plot
bigsources = repmat([""], length(supp.subsource), 1);
for i = 1:length(bigsources)
    if startsWith(supp.subsource{i},"real-")
        bigsources(i) = "Real data";
    elseif startsWith(supp.subsource{i},"reallike-")
        bigsources(i) = "Real-like.";
    elseif startsWith(supp.subsource{i},"manhat-")
        bigsources(i) = "Grid-based";
    elseif startsWith(supp.subsource{i},"random-")
        bigsources(i) = "Uniform random";
    elseif startsWith(supp.subsource{i},"recombined-")
        bigsources(i) = "Hybrid";
    elseif startsWith(supp.subsource{i}, "flowcluster-")
        bigsources(i) = "Flowcluster";
    else
        bigsources(i) = "Other instances";
    end
end
bigsourcescat = categorical(bigsources);

drawSources(model.pilot.Z, bigsourcescat, cmap);
title('Instance Categories')
print(gcf,'-dpng',[outputdir 'extone_sources.png']);
print(gcf,'-depsc',[outputdir 'extone_sources.eps']);

% flow cluster
fcsources = repmat([""], length(subs), 1);
for i = 1:length(fcsources)
    if startsWith(subs{i},"flowcluster-dhyper-fcycle")
        fcsources(i) = "Hcube x Triangle";
    elseif startsWith(subs{i},"flowcluster-dhyper-ftree")
        fcsources(i) = "Hcube x Tree";
    elseif startsWith(subs{i},"flowcluster-dhyper-fsquare")
        fcsources(i) = "Hcube x Square";
    elseif startsWith(subs{i},"flowcluster-ddrez-fcycle")
        fcsources(i) = "Drexx x Triangle";
    elseif startsWith(subs{i},"flowcluster-ddrez-ftree")
        fcsources(i) = "Drexx x Tree";
    elseif startsWith(subs{i},"flowcluster-ddrez-fsquare")
        fcsources(i) = "Drexx x Square";
    else
        fcsources(i) = "";
    end
end
fcsourcescat = categorical(fcsources);
typs = {"Hcube x Triangle", "Hcube x Tree", "Hcube x Square","Drexx x Triangle", "Drexx x Tree", "Drexx x Square"};
%typs = typs(1:2:5,2:4:6);

drawSources(model.pilot.Z, fcsourcescat, cmap, typs);
title('Flow-Cluster Instances')
print(gcf,'-dpng',[outputdir 'extone_flowcluster.png']);
print(gcf,'-depsc',[outputdir 'extone_flowcluster.eps']);

% hybrid
hybsources = repmat([""], length(subs), 1);
for i = 1:length(fcsources)
    if startsWith(subs{i},"recombined-ddrez")
        hybsources(i) = "DreXX";
    elseif startsWith(subs{i},"recombined-deucl")
        hybsources(i) = "Euclidean";
    elseif startsWith(subs{i},"recombined-dhypr")
        hybsources(i) = "Hypercube";
    elseif startsWith(subs{i},"recombined-dmanh")
        hybsources(i) = "Manhattan";
    elseif startsWith(subs{i},"recombined-dpalu")
        hybsources(i) = "Palubeckis";
    elseif startsWith(subs{i},"recombined-drand")
        hybsources(i) = "Random";
    elseif startsWith(subs{i},"recombined-dterm")
        hybsources(i) = "Terminal";
    else
        hybsources(i) = "";
    end
end
hybsourcescat = categorical(hybsources);
typs = {"DreXX", "Euclidean", "Hypercube", "Manhattan", "Palubeckis", "Random", "Terminal"};
%typs = typs(1:2:5,2:4:6);

drawSources(model.pilot.Z, hybsourcescat, cmap, typs);
title('Hybrid Instances')
print(gcf,'-dpng',[outputdir 'extone_hybrid.png']);
print(gcf,'-depsc',[outputdir 'extone_hybrid.eps']);

% specific plot
spec1sources = repmat([""], length(supp.subsource), 1);
for i = 1:length(spec1sources)
    if startsWith(supp.subsource{i},"other-gen-palubeckis")
        spec1sources(i) = "Palubeckis";
    elseif startsWith(supp.subsource{i},"terminal-gen")
        spec1sources(i) = "Terminal";
    elseif startsWith(supp.subsource{i},"hypercube")
        spec1sources(i) = "Hypercube";
    elseif startsWith(supp.subsource{i},"qapsat-gen")
        spec1sources(i) = "QAPSAT";
    elseif startsWith(supp.subsource{i},"other-drezner")
        spec1sources(i) = "DreXX";
    else
        spec1sources(i) = "";
    end
end
spec1sourcescat = categorical(spec1sources);
typs = {"Palubeckis", "Terminal", "Hypercube","QAPSAT", "DreXX"};

drawSources(model.pilot.Z, spec1sourcescat, cmap, typs);
title('Selected sub-classes')
print(gcf,'-dpng',[outputdir 'extone_specific.png']);
print(gcf,'-depsc',[outputdir 'extone_specific.eps']);


% reallike plot
RLsources = repmat([""], length(supp.subsource), 1);
for i = 1:length(RLsources)
    if startsWith(supp.subsource{i},"reallike-SF-euc-plu")
        RLsources(i) = "SFgen, StructPlus flows";
    elseif startsWith(supp.subsource{i},"reallike-SF-euc-ran")
        RLsources(i) = "SFgen, Random flows";
    elseif startsWith(supp.subsource{i},"reallike-SF-euc-str")
        RLsources(i) = "SFgen, Structured flows";
    elseif startsWith(supp.subsource{i},"reallike-gen-taiBN")
        RLsources(i) = "Tgen, normal distances";
    elseif startsWith(supp.subsource{i},"reallike-gen-taiBT")
        RLsources(i) = "Tgen, tilted distances";
    elseif startsWith(supp.subsource{i},"reallike-qaplib")
        RLsources(i) = "QAPLIB instances";
    else
        RLsources(i) = "";
    end
end
RLsourcescat = categorical(RLsources);
typs = {"SFgen, StructPlus flows", "SFgen, Random flows", "SFgen, Structured flows","Tgen, normal distances", "Tgen, tilted distances", "QAPLIB instances"};

drawSources(model.pilot.Z, RLsourcescat, cmap, typs);
title('Real-Like Instances')
print(gcf,'-dpng',[outputdir 'extone_reallike.png']);
print(gcf,'-depsc',[outputdir 'extone_reallike.eps']);

%manhattan plot
MHsources = repmat([""], length(supp.subsource), 1);
for i = 1:length(MHsources)
    if startsWith(supp.subsource{i},"manhat-gen-SF-plu")
        MHsources(i) = "SFgen, StructPlus flows";
    elseif startsWith(supp.subsource{i},"manhat-gen-SF-ran")
        MHsources(i) = "SFgen, Random flows";
    elseif startsWith(supp.subsource{i},"manhat-gen-SF-str")
        MHsources(i) = "SFgen, Structured flows";
    elseif startsWith(supp.subsource{i},"manhat-qaplib")
        MHsources(i) = "QAPLIB instances";
    else
        MHsources(i) = "";
    end
end
MHsourcescat = categorical(MHsources);
typs = {"SFgen, StructPlus flows", "SFgen, Random flows", "SFgen, Structured flows", "QAPLIB instances"};

drawSources(model.pilot.Z, MHsourcescat, cmap, typs);
title('Instances with Manhattan grid distances')
print(gcf,'-dpng',[outputdir 'extone_manhat.png']);
print(gcf,'-depsc',[outputdir 'extone_manhat.eps']);

%other1 plot
O1sources = repmat([""], length(supp.subsource), 1);
for i = 1:length(O1sources)
    if startsWith(supp.subsource{i},"terminal-gen")
        O1sources(i) = "Terminal";
    elseif startsWith(supp.subsource{i},"hypercube-gen")
        O1sources(i) = "Hypercube";
    elseif startsWith(supp.subsource{i},"other-gen-palubeckis")
        O1sources(i) = "Palubeckis";
    elseif startsWith(supp.subsource{i},"qapsat-gen-easy")
        O1sources(i) = "QAPSAT, easy";
    elseif startsWith(supp.subsource{i},"qapsat-gen-hard")
        O1sources(i) = "QAPSAT, hard";
    else
        O1sources(i) = "";
    end
end
O1sourcescat = categorical(O1sources);
typs = {"Terminal", "Hypercube", "Palubeckis", "QAPSAT, easy", "QAPSAT, hard"};

drawSources(model.pilot.Z, O1sourcescat, cmap, typs);
title('Other generated instances')
print(gcf,'-dpng',[outputdir 'extone_other1.png']);
print(gcf,'-depsc',[outputdir 'extone_other1.eps']);

%other2 plot
O2sources = repmat([""], length(supp.subsource), 1);
for i = 1:length(O2sources)
    if startsWith(supp.subsource{i},"other-drezner")
        O2sources(i) = "DreXX";
    elseif startsWith(supp.subsource{i},"other-qaplib-chr")
        O2sources(i) = "chr* (QAPLIB)";
    elseif startsWith(supp.subsource{i},"other-qaplib-lipa")
        O2sources(i) = "lipa* (QAPLIB)";
    elseif startsWith(supp.subsource{i},"other-qaplib-taic")
        O2sources(i) = "tai64c (QAPLIB)";
    else
        O2sources(i) = "";
    end
end
O2sourcescat = categorical(O2sources);
typs = {"DreXX", "chr* (QAPLIB)", "lipa* (QAPLIB)", "tai64c (QAPLIB)"};

drawSources(model.pilot.Z, O2sourcescat, cmap, typs);
title('Other benchmark instances')
print(gcf,'-dpng',[outputdir 'extone_other2.png']);
print(gcf,'-depsc',[outputdir 'extone_other2.eps']);

% algorithm performance plots
diffYraw = model.data.Yraw(:,2) - model.data.Yraw(:,1);

[h1, h2, h3] = drawScatterYraw(model.pilot.Z, diffYraw, "Actual Algorithm Performance", cmap);
legend([h1,h2,h3], ["BMA stronger", "Similar performance", "MMAS stronger"], 'Location', 'SouthOutside', 'NumColumns', 3);
print(gcf,'-dpng',[outputdir 'extone_realperf.png']);
print(gcf,'-depsc',[outputdir 'extone_realperf.eps']);

% SVM plots
hold off
clf
drawBinaryPerformance(model.pilot.Z, model.pythia.Yhat(:,1), "SVM prediction of BMA performance", cmap)
print(gcf,'-dpng',[outputdir 'extone_svm_bma.png']);
print(gcf,'-depsc',[outputdir 'extone_svm_bma.eps']);
clf
drawBinaryPerformance(model.pilot.Z, model.pythia.Yhat(:,2), "SVM prediction of MMAS performance", cmap)
print(gcf,'-dpng',[outputdir 'extone_svm_mmas.png']);
print(gcf,'-depsc',[outputdir 'extone_svm_mmas.eps']);

% feature plots
Xaux = (model.data.X-min(model.data.X,[],1))./range(model.data.X,1);
longfeat = {'Distance Normalised Mean', "Flow Normalised Mean", "Distance Sparsity", "Flow Dominance", "Distance Skewness", "Gilmore Lawler Bound", "Escape Probability"};
for i=1:nfeats
    clf;
    drawScatter(model.pilot.Z, Xaux(:,i),...
                strrep(model.data.featlabels{i},'_',' '), cmap);
    % line(model.cloist.Zedge(:,1), model.cloist.Zedge(:,2), 'LineStyle', '-', 'Color', 'r');
    title(longfeat{i});
    print(gcf,'-dpng',[outputdir 'extone_feature_' model.data.featlabels{i} '.png']);
    print(gcf,'-depsc',[outputdir 'extone_feature_' model.data.featlabels{i} '.eps']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if false

%figure
clf
%set(gcf, 'Position',  [0, 100, 800, 800])
%drawSubSources(model.pilot.Z, subS, model.data.S, rootdir);
print(gcf,'-depsc',[rootdir 'subsource.eps']);


clf
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );
A = model.pilot.A;
feats = size(model.pilot.A,2);
[X,I] = sort(atan2(A(1,:),A(2,:)));
A = A(:,I);
A = A(:,[1:2:end,2:2:end]);
A = A(:,end:-1:1);
bnd = max(abs(A(:))) * 1.2;
hold on
for i = 1:feats
    x = A(:,i);
    c = hsv2rgb([i/feats, 1, 0.6]);
    drawArrow([0,x(1)],[0,x(2)],'linewidth',3,'color',c);
end
axis square;
xlim([-bnd bnd]);
ylim([-bnd bnd]);
%set(gcf, 'Position',  [0, 100, 800, 800])
hold off
print(gcf,'-dpng',[rootdir 'arrows.png']);
print(gcf,'-depsc',[rootdir 'arrows.eps']);

clf
Z = model.pilot.Z;
cats = model.data.Ybin*[1;2];
BMAonly = (cats == 1);
MMASonly = (cats == 2);
allalgs = (cats == 3);
hold on
scatter(Z(allalgs,1), Z(allalgs,2),10,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
scatter(Z(MMASonly,1), Z(MMASonly,2),80,'Marker','pentagram','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 1], 'LineWidth', 2);
scatter(Z(BMAonly,1), Z(BMAonly,2),80,'Marker','pentagram','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0], 'LineWidth', 2);
hold off
xlabel('z_{1}'); ylabel('z_{2}');
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([floor(min(Z(:,1))) ceil(max(Z(:,1))) floor(min(Z(:,2))) ceil(max(Z(:,2)))]);
legend("All similar","Good MMAS only","Good BMA only","Location","eastoutside")
print(gcf,'-depsc',[rootdir 'distribution_portfolio_alt.eps']);
print(gcf,'-dpng',[rootdir 'distribution_portfolio_alt.png']);

lineqns = zeros(2,size(model.data.Yraw,2));
for i = 1:size(model.data.Yraw,2)
    clf
    Ydata = model.data.Yraw(:,i);
    Pr0hat = model.pythia.Pr0hat(:,i);
    X = [ones(size(Pr0hat,1),1),Pr0hat];
    Leq = X \ Ydata;
    %pred = Leq(1) + Leq(2)*Yraw;
    scatter(Pr0hat, Ydata)
    hold on
    fplot(@(x) Leq(1) + Leq(2)*x, [min(Pr0hat),max(Pr0hat)]);
    xlabel("Pr0hat from SVM")
    ylabel(strcat("Performance of ", model.data.algolabels{i}));
    title("Linear Regression Relation between Pr0hat and Performance")
    hold off
    print(gcf,'-dpng',[rootdir model.data.algolabels{i} '_linreg.png']);

    lineqns(1,i) = Leq(1);
    lineqns(2,i) = Leq(2);
end

clf
proj_perform = model.pythia.Pr0hat .* lineqns(2,:) + lineqns(1,:);
[~, idx] = min(proj_perform,[],2);
hold on
Z = model.pilot.Z;
for i = 1:size(model.data.Yraw,2)
    choice(:,i) = (idx == i);
    color = [0 0 0]; % this is a hack but im in a hurry
    color(i) = 1;
    color = color * [1 0 0; 0 0 1; 0 1 0];
    scatter(Z(choice(:,i),1), Z(choice(:,i),2),20,'MarkerEdgeColor',color,'MarkerFaceColor',color);
end
title("Pr0hat informed prediction of best alg");
legend("BMA","MMAS");
print(gcf,'-dpng',[rootdir 'Pr0hat_prediction.png']);
hold off

end


function handle = drawSubSources(Z, subS, supS, rootdir)

clf
ubound = ceil(max(Z));
lbound = floor(min(Z));
sourcelabels = cellstr(unique(subS));
nsources = length(sourcelabels);
clrs = flipud(colorcube(nsources+4));
clrs = clrs(2:end-3,:);
clrs = clrs([1:4:nsources,2:4:nsources,3:4:nsources,4:4:nsources],:);
handle = zeros(nsources,1);

supercats = categories(supS);
%symbols = {'o','x','s', '^', '+', "pentagram", '<'};
symbols = {'>','+','s','*', '^', 'o', "pentagram", 'hexagram', 'd'};
sizes = [5,5,5,5,5,5,5,5,5];

for i=nsources:-1:1
    subsI = (subS==sourcelabels{i});
    super = supS(subsI);
    superI = find(super(1) == supercats);
    
    handle(i) = line(Z(subsI,1), ...
                     Z(subsI,2), ...
                     'LineStyle', 'none', ...
                     'Marker', symbols{superI}, ...
                     'Color', clrs(i,:), ...
                     'MarkerFaceColor', clrs(i,:), ...
                     'MarkerSize', sizes(superI));
end
xlabel('z_{1}'); ylabel('z_{2}'); title('Sub-Sources');
legend(handle, sourcelabels, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

print(gcf,'-dpng',[rootdir 'distribution_subsources.png']);

end