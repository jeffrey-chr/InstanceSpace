function QAP_makeEPS(container,rootdir,outdir)

% -------------------------------------------------------------------------
% pgnscript.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------
% Modified version to produce .EPS files
% -------------------------------------------------------------------------
% Preliminaries
scriptfcn;
colormap('parula');
nfeats = size(container.data.X,2);
nalgos = size(container.data.Y,2);
Xaux = (container.data.X-min(container.data.X,[],1))./range(container.data.X,1);
Yind = (container.data.Yraw-min(container.data.Yraw,[],1))./range(container.data.Yraw,1);
Yglb = log10(container.data.Yraw+1);
Yglb = (Yglb-min(Yglb(:)))./range(Yglb(:));
if container.opts.trace.usesim
    Yfoot = container.pythia.Yhat;
    Pfoot = container.pythia.selection0;
else
    Yfoot = container.data.Ybin;
    Pfoot = container.data.P;
end

xl = [-4 4];
yl = [-4 4];
% -------------------------------------------------------------------------
disp('=========================================================================');
disp('-> Producing the plots.');
% -------------------------------------------------------------------------
% Drawing feature plots
for i=1:nfeats
    clf;
    drawScatter(container.pilot.Z, Xaux(:,i),...
                strrep(container.data.featlabels{i},'_',' '));
    % line(model.cloist.Zedge(:,1), model.cloist.Zedge(:,2), 'LineStyle', '-', 'Color', 'r');
    title(''); xlim(xl); ylim(yl);
    set(findall(gcf,'-property','SizeData'),'SizeData',60);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    print(gcf,'-depsc',[outdir 'eps\distribution_feature_' container.data.featlabels{i} '.eps']);
    print(gcf,'-dpng',[outdir 'eps\distribution_feature_' container.data.featlabels{i} '.png']);
    %print([outdir 'distribution_feature_' container.data.featlabels{i} '.eps'], '-depsc');
end
% -------------------------------------------------------------------------
% Drawing algorithm performance/footprint plots
for i=1:nalgos
    % Actual performance, normalized globaly
    clf;
    drawScatter(container.pilot.Z, Yglb(:,i), ...
                strrep(container.data.algolabels{i},'_',' '));
    print(gcf,'-depsc',[outdir 'eps\distribution_performance_global_normalized_' container.data.algolabels{i} '.eps']);
    % Actual performance, normalized individualy
    clf;
    drawScatter(container.pilot.Z, Yind(:,i), ...
                strrep(container.data.algolabels{i},'_',' '));
    print(gcf,'-depsc',[outdir 'eps\distribution_performance_individual_normalized_' container.data.algolabels{i} '.eps']);
    % Actual binary performance
    clf;
    drawBinaryPerformance(container.pilot.Z, container.data.Ybin(:,i), ...
                          strrep(container.data.algolabels{i},'_',' '));
    title(''); xlim(xl); ylim(yl);
    set(findall(gcf,'-property','MarkerSize'),'MarkerSize',15);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    s=findobj('type','legend');
    delete(s);
    saveas(gcf, [outdir 'eps\binary_performance_' container.data.algolabels{i} '.fig'])
    print(gcf,'-depsc',[outdir 'eps\binary_performance_' container.data.algolabels{i} '.eps']);
    print(gcf,'-dpng',[outdir 'eps\binary_performance_' container.data.algolabels{i} '.png']);
    % Drawing the SVM's predictions of good performance
    try
        clf;
        drawBinaryPerformance(container.pilot.Z, container.pythia.Yhat(:,i), ...
                              strrep(container.data.algolabels{i},'_',' '));
        print(gcf,'-depsc',[outdir 'eps\binary_svm_' container.data.algolabels{i} '.eps']);
    catch
        disp('No SVM model has been trained');
    end
    % Drawing the footprints for good and bad performance acording to the
    % binary measure 
    try 
        clf;
        drawGoodBadFootprint(container.pilot.Z, ...
                             container.trace.good{i}, ...
                             Yfoot(:,i), ...
                             strrep(container.data.algolabels{i},'_',' '));
        print(gcf,'-depsc',[outdir 'eps\footprint_' container.data.algolabels{i} '.eps']);
    catch
        disp('No Footprint has been calculated');
    end
end
% ---------------------------------------------------------------------
% Plotting the number of good algos
clf;
drawScatter(container.pilot.Z, container.data.numGoodAlgos./nalgos, 'Percentage of good algorithms');
print(gcf,'-depsc',[outdir 'eps\distribution_number_good_algos.eps']);
% ---------------------------------------------------------------------
% Drawing the algorithm performance
clf;
drawPortfolioSelections(container.pilot.Z, container.data.P, container.data.algolabels, 'Best algorithm');
print(gcf,'-depsc',[outdir 'eps\distribution_portfolio.eps']);
% ---------------------------------------------------------------------
% Drawing the SVM's recommendations
clf;
drawPortfolioSelections(container.pilot.Z, container.pythia.selection0, container.data.algolabels, 'Predicted best algorithm');
print(gcf,'-depsc',[outdir 'eps\distribution_svm_portfolio.eps']);
% ---------------------------------------------------------------------
% Drawing the footprints as portfolio.
clf;
drawPortfolioFootprint(container.pilot.Z, container.trace.best, Pfoot, container.data.algolabels);
title(''); xlim(xl); ylim(yl);
s=findobj('type','legend');
delete(s);
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',15);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(gcf, [outdir 'eps\footprint_portfolio.fig'])
print(gcf,'-depsc',[outdir 'eps\footprint_portfolio.eps']);
print(gcf,'-dpng',[outdir 'eps\footprint_portfolio.png']);
% ---------------------------------------------------------------------
% Plotting the model.data.beta score
clf;
drawBinaryPerformance(container.pilot.Z, container.data.beta, '\beta score');
print(gcf,'-depsc',[outdir 'eps\distribution_beta_score.eps']);
% ---------------------------------------------------------------------
% Drawing the sources of the instances if available
if isfield(container.data,'S')
    clf;
    drawSources(container.pilot.Z, container.data.S);
    print(gcf,'-depsc',[outdir 'eps\distribution_sources.eps']);
end

% draw which instances are in QAPLIB
clf;
libreg = "^rou.*$|^tai.*a.*$|^nug.*$|^sko.*$|^wil.*$|^had.*$|^scr.*$|^tho.*$|^ste.*$|^els.*$|^bur.*$|^chr.*$|^kra.*$|^esc.*$|^tai.*b.*$|^lipa.*$|^tai.*c.*$";
libsource = [];
for in = 1:length(container.data.instlabels)
    aonly = char(container.data.instlabels(in));
    if ~isempty(regexp(aonly,libreg,'ONCE'))
        libsource = [libsource;"QAPLIB"];
    else
        libsource = [libsource;"OTHER"];
    end
end
libsource = categorical(libsource);

if true
    Z = container.pilot.Z;
    S = libsource;
    ubound = ceil(max(Z));
    lbound = floor(min(Z));
    sourcelabels = cellstr(unique(S));
    nsources = length(sourcelabels);
    %clrs = flipud(lines(nsources));
    clrs = [0.6 0.6 0.6; 1 0 0];
    sizes = [15 15];
    markers = ['.' 'p'];
    handle = zeros(nsources,1);
    for i=nsources:-1:1
        handle(i) = line(Z(S==sourcelabels{i},1), ...
                         Z(S==sourcelabels{i},2), ...
                         'LineStyle', 'none', ...
                         'Marker', markers(i), ...
                         'Color', clrs(i,:), ...
                         'MarkerFaceColor', clrs(i,:), ...
                         'MarkerSize', sizes(i));
    end
    xlabel('z_{1}'); ylabel('z_{2}'); title('Sources');
    legend(handle, sourcelabels, 'Location', 'NorthEastOutside');
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
    axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);
end

title(''); xlim(xl); ylim(yl);
s=findobj('type','legend');
delete(s);
%set(findall(gcf,'-property','MarkerSize'),'MarkerSize',15);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
print(gcf,'-depsc',[outdir 'eps\distribution_qaplib.eps']);
print(gcf,'-dpng',[outdir 'eps\distribution_qaplib.png']);

% draw which instances are in QAPLIB
clf;
manreg = "^nug.*$|^sko.*$|^wil.*$|^had.*$|^scr.*$|^tho.*$|^stf.*m[prs].*$";
realreg = "^ste.*$|^els.*$|^bur.*$|^kra.*$|^esc.*$";
likereg = "^tai.*b.*$|^stf.*e[prs].*|^xtab.*$";
typesource = [];
for in = 1:length(container.data.instlabels)
    aonly = char(container.data.instlabels(in));
    if ~isempty(regexp(aonly,manreg,'ONCE'))
        typesource = [typesource;"MANHATTAN"];
    elseif ~isempty(regexp(aonly,realreg,'ONCE'))
        typesource = [typesource;"REAL"];
    elseif ~isempty(regexp(aonly,likereg,'ONCE'))
        typesource = [typesource;"REALLIKE"];
    else
        typesource = [typesource;"ZZZ"];
    end
end
typesource = categorical(typesource);

if true
    Z = container.pilot.Z;
    S = typesource;
    ubound = ceil(max(Z));
    lbound = floor(min(Z));
    sourcelabels = cellstr(unique(S));
    nsources = length(sourcelabels);
    %clrs = flipud(lines(nsources));
    clrs = [1 0 0; 0 1 0; 0 0 1; 0.6 0.6 0.6];
    sizes = [15 15 15 15];
    markers = ['s' 'p' '^' '.'];
    handle = zeros(nsources,1);
    for i=nsources:-1:1
        handle(i) = line(Z(S==sourcelabels{i},1), ...
                         Z(S==sourcelabels{i},2), ...
                         'LineStyle', 'none', ...
                         'Marker', markers(i), ...
                         'Color', clrs(i,:), ...
                         'MarkerFaceColor', clrs(i,:), ...
                         'MarkerSize', sizes(i));
    end
    xlabel('z_{1}'); ylabel('z_{2}'); title('Sources');
    legend(handle, sourcelabels, 'Location', 'NorthEastOutside');
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
    axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);
end

title(''); xlim(xl); ylim(yl);
s=findobj('type','legend');
delete(s);
%set(findall(gcf,'-property','MarkerSize'),'MarkerSize',15);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
print(gcf,'-depsc',[outdir 'eps\distribution_type.eps']);
print(gcf,'-dpng',[outdir 'eps\distribution_type.png']);

clf;
if true
    Z = container.pilot.Z;
    S = typesource;
    ubound = ceil(max(Z));
    lbound = floor(min(Z));
    %clrs = flipud(lines(nsources));
    clrs = [1 0 0; 0 0 1; 0.6 0.6 0.6];
    sizes = [15 15 15];
    markers = ['x' '+' '.'];
    handle = zeros(3,1);
    cats = container.data.Ybin*[1;2];
    BMAonly = (cats == 1);
    MMASonly = (cats == 2);
    BMAandMMAS = (cats == 3);
    handle(3) = line(Z(BMAandMMAS,1), ...
                     Z(BMAandMMAS,2), ...
                     'LineStyle', 'none', ...
                     'Marker', markers(3), ...
                     'Color', clrs(3,:), ...
                     'MarkerFaceColor', clrs(3,:), ...
                     'MarkerSize', sizes(3));
    handle(1) = line(Z(BMAonly,1), ...
                     Z(BMAonly,2), ...
                     'LineStyle', 'none', ...
                     'Marker', markers(1), ...
                     'Color', clrs(1,:), ...
                     'MarkerFaceColor', clrs(1,:), ...
                     'MarkerSize', sizes(1));
    handle(2) = line(Z(MMASonly,1), ...
                     Z(MMASonly,2), ...
                     'LineStyle', 'none', ...
                     'Marker', markers(2), ...
                     'Color', clrs(2,:), ...
                     'MarkerFaceColor', clrs(2,:), ...
                     'MarkerSize', sizes(2));
    
    xlabel('z_{1}'); ylabel('z_{2}'); %title('Sources');
    %legend(handle, sourcelabels, 'Location', 'NorthEastOutside');
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
    axis square; axis([-4 4 -4 4]);
    print(gcf,'-depsc',[outdir 'eps\findbad.eps']);
end

lineqns = zeros(2,size(container.data.Yraw,2));
for i = 1:size(container.data.Yraw,2)
    clf
    Ydata = container.data.Yraw(:,i);
    Pr0hat = container.pythia.Pr0hat(:,i);
    X = [ones(size(Pr0hat,1),1),Pr0hat];
    Leq = X \ Ydata;
    %pred = Leq(1) + Leq(2)*Yraw;

    lineqns(1,i) = Leq(1);
    lineqns(2,i) = Leq(2);
end

clf
proj_perform = container.pythia.Pr0hat .* lineqns(2,:) + lineqns(1,:);
[~, idx] = min(proj_perform,[],2);
hold on
Z = container.pilot.Z;
i=1;
choice(:,i) = (idx == i);
    color = [1 0 0]; % this is a hack but im in a hurry
    scatter(Z(choice(:,i),1), Z(choice(:,i),2),20,'MarkerEdgeColor',color,'MarkerFaceColor',color);
i=2;
   choice(:,i) = (idx == i);
    color = [0 0 1]; % this is a hack but im in a hurry
    scatter(Z(choice(:,i),1), Z(choice(:,i),2),20,'MarkerEdgeColor',color,'MarkerFaceColor',color);

%title("Pr0hat informed prediction of best alg");
%legend("BMA","MMAS");
print(gcf,'-depsc',[outdir 'eps\Pr0hat_prediction.eps']);
hold off