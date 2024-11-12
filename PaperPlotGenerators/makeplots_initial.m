% Create plots of intitial space 
clear;

rootdir = '..\QAPdata\';
model = load([rootdir 'model.mat']);
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});

outputdir = '.\output\';

nfeats = length(model.data.featlabels);

f = gcf;
goodpos = [50 50 800 600];
f.Position = goodpos;

cmap = @copper;

alteredscriptfcn

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
print(gcf,'-dpng',[outputdir 'init_qaplib.png']);
print(gcf,'-depsc',[outputdir 'init_qaplib.eps']);

% category plot
bigsources = repmat([""], length(supp.subsource), 1);
for i = 1:length(bigsources)
    if startsWith(supp.subsource{i},"real-")
        bigsources(i) = "Real data";
    elseif startsWith(supp.subsource{i},"reallike-")
        bigsources(i) = "Real-like";
    elseif startsWith(supp.subsource{i},"manhat-")
        bigsources(i) = "Grid-based";
    elseif startsWith(supp.subsource{i},"random-")
        bigsources(i) = "Uniform random";
    else
        bigsources(i) = "Other";
    end
end
bigsourcescat = categorical(bigsources);

drawSources(model.pilot.Z, bigsourcescat, cmap);
title('Instance Categories')
print(gcf,'-dpng',[outputdir 'init_sources.png']);
print(gcf,'-depsc',[outputdir 'init_sources.eps']);

% reallike plot
RLsources = repmat([""], length(supp.subsource), 1);
for i = 1:length(RLsources)
    if startsWith(supp.subsource{i},"reallike-SF-euc-plu")
        RLsources(i) = "SFgen, Struct+";
    elseif startsWith(supp.subsource{i},"reallike-SF-euc-ran")
        RLsources(i) = "SFgen, Random";
    elseif startsWith(supp.subsource{i},"reallike-SF-euc-str")
        RLsources(i) = "SFgen, Struct";
    elseif startsWith(supp.subsource{i},"reallike-gen-taiBN")
        RLsources(i) = "Taigen, Normal";
    elseif startsWith(supp.subsource{i},"reallike-gen-taiBT")
        RLsources(i) = "Taigen, Tilt";
    elseif startsWith(supp.subsource{i},"reallike-qaplib")
        RLsources(i) = "tai*b (QAPLIB)";
    else
        RLsources(i) = "";
    end
end
RLsourcescat = categorical(RLsources);
typs = {"SFgen, Struct+", "SFgen, Random", "SFgen, Struct","Tgen, Normal", "Tgen, Tilt", "tai*b (QAPLIB)"};

drawSources(model.pilot.Z, RLsourcescat, cmap, typs);
title('Real-Like Instances')
print(gcf,'-dpng',[outputdir 'init_reallike.png']);
print(gcf,'-depsc',[outputdir 'init_reallike.eps']);

%manhattan plot
MHsources = repmat([""], length(supp.subsource), 1);
for i = 1:length(MHsources)
    if startsWith(supp.subsource{i},"manhat-gen-SF-plu")
        MHsources(i) = "SFgen, Struct+";
    elseif startsWith(supp.subsource{i},"manhat-gen-SF-ran")
        MHsources(i) = "SFgen, Random";
    elseif startsWith(supp.subsource{i},"manhat-gen-SF-str")
        MHsources(i) = "SFgen, Struct";
    elseif startsWith(supp.subsource{i},"manhat-qaplib")
        MHsources(i) = "QAPLIB";
    else
        MHsources(i) = "";
    end
end
MHsourcescat = categorical(MHsources);
typs = {"SFgen, Struct+", "SFgen, Random", "SFgen, Struct", "QAPLIB"};

drawSources(model.pilot.Z, MHsourcescat, cmap, typs);
title('Instances with Manhattan grid distances')
print(gcf,'-dpng',[outputdir 'init_manhat.png']);
print(gcf,'-depsc',[outputdir 'init_manhat.eps']);

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
print(gcf,'-dpng',[outputdir 'init_other1.png']);
print(gcf,'-depsc',[outputdir 'init_other1.eps']);

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
print(gcf,'-dpng',[outputdir 'init_other2.png']);
print(gcf,'-depsc',[outputdir 'init_other2.eps']);

% algorithm performance plots
diffYraw = model.data.Yraw(:,2) - model.data.Yraw(:,1);

[h1, h2, h3] = drawScatterYraw(model.pilot.Z, diffYraw, "Actual Algorithm Performance", cmap);
perfleg = legend([h1,h2,h3], ["BMA stronger", "Similar performance", "MMAS stronger"], 'Location', 'SouthOutside', 'NumColumns', 3);
cpos = perfleg.Position;
perfleg.set("Position", [(1 - cpos(3))/2, 0.01, cpos(3), cpos(4)]);
gca().set("Position", gca().Position + [0.02 0.04 -0.04 -0.04]);
print(gcf,'-dpng',[outputdir 'init_realperf.png']);
print(gcf,'-depsc',[outputdir 'init_realperf.eps']);

close all
f = gcf;
f.Position = goodpos;

% feature plots
Xaux = (model.data.X-min(model.data.X,[],1))./range(model.data.X,1);
longfeat = {'Distance Sparsity', "Distance Triangle Ineq. Sat.", "Distance Betafit Alpha", "Distance Near Similarity", "Cumulative Integral", "Average Distance to Optima"};
for i=1:nfeats
    clf;
    drawScatter(model.pilot.Z, Xaux(:,i),...
                strrep(model.data.featlabels{i},'_',' '), cmap);
    % line(model.cloist.Zedge(:,1), model.cloist.Zedge(:,2), 'LineStyle', '-', 'Color', 'r');
    title(longfeat(i));
    print(gcf,'-dpng',[outputdir 'init_feature_' model.data.featlabels{i} '.png']);
    print(gcf,'-depsc',[outputdir 'init_feature_' model.data.featlabels{i} '.eps']);
end

% SVM plots
hold off
clf
drawBinaryPerformance(model.pilot.Z, model.pythia.Yhat(:,1), "SVM prediction of BMA performance", cmap)
print(gcf,'-dpng',[outputdir 'init_svm_bma.png']);
print(gcf,'-depsc',[outputdir 'init_svm_bma.eps']);
clf
drawBinaryPerformance(model.pilot.Z, model.pythia.Yhat(:,2), "SVM prediction of MMAS performance", cmap)
print(gcf,'-dpng',[outputdir 'init_svm_mmas.png']);
print(gcf,'-depsc',[outputdir 'init_svm_mmas.eps']);

%default SVM prediction
clf
drawPortfolioSelections(model.pilot.Z, model.pythia.selection0, model.data.algolabels, "Combined SVM prediction based on global precision", cmap)
print(gcf,'-dpng',[outputdir 'init_svm_combined1.png']);
print(gcf,'-depsc',[outputdir 'init_svm_combined1.eps']);

basicgood = 0;
basicbest = 0;
for i = 1:size(model.pythia.selection0, 1)
    s = model.pythia.selection0(i);
    if (s > 0)
        if (model.data.Ybin(i,s) == 1)
            basicgood = basicgood + 1;
        end
        chosen = model.data.Yraw(i,s);
        notchosen = model.data.Yraw(i,3-s);
        if (chosen <= notchosen)
            basicbest = basicbest+1;
        end
    end
end

%improved SVM prediction
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
% hold on
% Z = model.pilot.Z;
% colors = cmap(2);
% for i = 1:size(model.data.Yraw,2)
%     choice(:,i) = (idx == i);
%     color = [0 0 0]; % this is a hack but im in a hurry
%     color(i) = 1;
%     %color = color * [1 0 0; 0 0 1]';
%     scatter(Z(choice(:,i),1), Z(choice(:,i),2),20,'MarkerEdgeColor',color,'MarkerFaceColor',color);
% end

newgood = 0;
newbest = 0;
for i = 1:size(model.pythia.selection0, 1)
    s = idx(i);
    if (s > 0)
        if (model.data.Ybin(i,s) == 1)
            newgood = newgood + 1;
        end
        chosen = model.data.Yraw(i,s);
        notchosen = model.data.Yraw(i,3-s);
        if (chosen <= notchosen)
            newbest = newbest+1;
        end
    end
end
basicgood = basicgood / size(model.pythia.selection0, 1);
basicbest = basicbest / size(model.pythia.selection0, 1);
newgood = newgood / size(model.pythia.selection0, 1);
newbest = newbest / size(model.pythia.selection0, 1);

drawPortfolioSelections(model.pilot.Z, idx, model.data.algolabels, "Combined SVM prediction based on local confidence", cmap)
title("Combined SVM prediction based on local confidence");
%legend("BMA","MMAS");
print(gcf,'-dpng',[outputdir 'init_svm_combined2.png']);
print(gcf,'-depsc',[outputdir 'init_svm_combined2.eps']);
hold off

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