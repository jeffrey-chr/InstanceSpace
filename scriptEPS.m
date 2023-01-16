function scriptEPS(container,rootdir)
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
    print(gcf,'-depsc',[rootdir 'eps\distribution_feature_' container.data.featlabels{i} '.eps']);
    print(gcf,'-dpng',[rootdir 'eps\distribution_feature_' container.data.featlabels{i} '.png']);
end
% -------------------------------------------------------------------------
% Drawing algorithm performance/footprint plots
for i=1:nalgos
    % Actual performance, normalized globaly
    clf;
    drawScatter(container.pilot.Z, Yglb(:,i), ...
                strrep(container.data.algolabels{i},'_',' '));
    print(gcf,'-depsc',[rootdir 'eps\distribution_performance_global_normalized_' container.data.algolabels{i} '.eps']);
    % Actual performance, normalized individualy
    clf;
    drawScatter(container.pilot.Z, Yind(:,i), ...
                strrep(container.data.algolabels{i},'_',' '));
    print(gcf,'-depsc',[rootdir 'eps\distribution_performance_individual_normalized_' container.data.algolabels{i} '.eps']);
    % Actual binary performance
    clf;
    drawBinaryPerformance(container.pilot.Z, container.data.Ybin(:,i), ...
                          strrep(container.data.algolabels{i},'_',' '));
    title(''); xlim(xl); ylim(yl);
    set(findall(gcf,'-property','MarkerSize'),'MarkerSize',15);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    s=findobj('type','legend');
    delete(s);
    saveas(gcf, [rootdir 'eps\binary_performance_' container.data.algolabels{i} '.fig'])
    print(gcf,'-depsc',[rootdir 'eps\binary_performance_' container.data.algolabels{i} '.eps']);
    print(gcf,'-dpng',[rootdir 'eps\binary_performance_' container.data.algolabels{i} '.png']);
    % Drawing the SVM's predictions of good performance
    try
        clf;
        drawBinaryPerformance(container.pilot.Z, container.pythia.Yhat(:,i), ...
                              strrep(container.data.algolabels{i},'_',' '));
        print(gcf,'-depsc',[rootdir 'eps\binary_svm_' container.data.algolabels{i} '.eps']);
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
        print(gcf,'-depsc',[rootdir 'eps\footprint_' container.data.algolabels{i} '.eps']);
    catch
        disp('No Footprint has been calculated');
    end
end
% ---------------------------------------------------------------------
% Plotting the number of good algos
clf;
drawScatter(container.pilot.Z, container.data.numGoodAlgos./nalgos, 'Percentage of good algorithms');
print(gcf,'-depsc',[rootdir 'eps\distribution_number_good_algos.eps']);
% ---------------------------------------------------------------------
% Drawing the algorithm performance
clf;
drawPortfolioSelections(container.pilot.Z, container.data.P, container.data.algolabels, 'Best algorithm');
print(gcf,'-depsc',[rootdir 'eps\distribution_portfolio.eps']);
% ---------------------------------------------------------------------
% Drawing the SVM's recommendations
clf;
drawPortfolioSelections(container.pilot.Z, container.pythia.selection0, container.data.algolabels, 'Predicted best algorithm');
print(gcf,'-depsc',[rootdir 'eps\distribution_svm_portfolio.eps']);
% ---------------------------------------------------------------------
% Drawing the footprints as portfolio.
clf;
drawPortfolioFootprint(container.pilot.Z, container.trace.best, Pfoot, container.data.algolabels);
title(''); xlim(xl); ylim(yl);
s=findobj('type','legend');
delete(s);
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',15);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
saveas(gcf, [rootdir 'eps\footprint_portfolio.fig'])
print(gcf,'-depsc',[rootdir 'eps\footprint_portfolio.eps']);
print(gcf,'-dpng',[rootdir 'eps\footprint_portfolio.png']);
% ---------------------------------------------------------------------
% Plotting the model.data.beta score
clf;
drawBinaryPerformance(container.pilot.Z, container.data.beta, '\beta score');
print(gcf,'-depsc',[rootdir 'eps\distribution_beta_score.eps']);
% ---------------------------------------------------------------------
% Drawing the sources of the instances if available
if isfield(container.data,'S')
    clf;
    drawSources(container.pilot.Z, container.data.S);
    print(gcf,'-depsc',[rootdir 'eps\distribution_sources.eps']);
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
print(gcf,'-depsc',[rootdir 'eps\distribution_qaplib.eps']);
print(gcf,'-dpng',[rootdir 'eps\distribution_qaplib.png']);

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
print(gcf,'-depsc',[rootdir 'eps\distribution_type.eps']);
print(gcf,'-dpng',[rootdir 'eps\distribution_type.png']);