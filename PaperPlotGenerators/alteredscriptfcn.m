function alteredscriptfcn
% -------------------------------------------------------------------------
% scriptfcn.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

writeArray2CSV = @(data,colnames,rownames,filename) writetable(array2table(data,'VariableNames',colnames,...
                                                                                'RowNames',rownames),...
                                                               filename,'WriteRowNames',true);
writeCell2CSV = @(data,colnames,rownames,filename) writetable(cell2table(data,'VariableNames',colnames,...
                                                                              'RowNames',rownames),...
                                                              filename,'WriteRowNames',true);
makeBndLabels = @(data) arrayfun(@(x) strcat('bnd_pnt_',num2str(x)),1:size(data,1),'UniformOutput',false);
colorscale  = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data,[],1)), range(data)));
colorscaleg = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data(:))), range(data(:))));

assignin('caller','writeArray2CSV',writeArray2CSV);
assignin('caller','writeCell2CSV',writeCell2CSV);
assignin('caller','makeBndLabels',makeBndLabels);
assignin('caller','colorscale',colorscale);
assignin('caller','colorscaleg',colorscaleg);
assignin('caller','drawSources',@drawSources);
assignin('caller','drawScatter',@drawScatter);
assignin('caller','drawScatterInfer',@drawScatterInfer);
assignin('caller','drawScatterYraw',@drawScatterYraw);
assignin('caller','drawPortfolioSelections',@drawPortfolioSelections);
assignin('caller','drawPortfolioFootprint',@drawPortfolioFootprint);
assignin('caller','drawGoodBadFootprint',@drawGoodBadFootprint);
assignin('caller','drawFootprint',@drawFootprint);
assignin('caller','drawBinaryPerformance',@drawBinaryPerformance);

end
% =========================================================================
% SUBFUNCTIONS
% =========================================================================
function handle = drawSources(Z, S, cmap, sourcelabels)

ubound = ceil(max(Z));
lbound = floor(min(Z));
if nargin < 4
    sourcelabels = cellstr(unique(S));
end
nsources = length(sourcelabels);
clrs = flipud(cmap(nsources));
clrs = [clrs(1:2:nsources,:); clrs(2:2:nsources,:)];
handle = zeros(nsources,1);
markertypes = [repmat('o', 1, length(1:2:nsources)), repmat('x', 1, length(2:2:nsources))];
markersizes = [repmat(20, 1, length(1:2:nsources)), repmat(20, 1, length(2:2:nsources))];
handle2 = zeros(nsources, 1);
handle3 = scatter(Z(isundefined(S),1), Z(isundefined(S),2), 5, [0.83 0.83 0.83], 'o', 'filled');
hold on
for i=nsources:-1:1
    % line(Z(S==sourcelabels{i},1), ...
    %      Z(S==sourcelabels{i},2), ...
    %      'LineStyle', 'none', ...
    %      'Marker', markertypes(rem(i-1,2)+1), ...
    %      'Color', clrs(i,:), ...
    %      'MarkerFaceColor', clrs(i,:), ...
    %      'MarkerSize', markersizes(rem(i-1,2)+1));
    
    handle2(i) = scatter(Z(S==sourcelabels{i},1), Z(S==sourcelabels{i},2), markersizes(i), clrs(i,:), markertypes(i));
    handle(i) = patch([0 0],[0 0], clrs(i,:), 'EdgeColor','none');
end
hold off
xlabel('z_{1}'); ylabel('z_{2}'); title('Sources');
leg = legend(handle2, sourcelabels, 'Location', 'SouthOutside', 'NumColumns', min(ceil(nsources/3),2));
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-0.5 ubound(1)+0.5 lbound(2)-0.5 ubound(2)+0.5]);
% cpos = leg.Position;
% leg.set("Position", [(1 - cpos(3))/2, 0.01, cpos(3), cpos(4)]);
% gca().set("Position", gca().Position + [0.03 0.06 -0.06 -0.06]);

end
% =========================================================================
function handle = drawScatterInfer(Z, X, titlelabel, cmap)

ubound = ceil(max(Z));
lbound = floor(min(Z));
handle = scatter(Z(:,1), Z(:,2), 8, X, 'filled');
colormap(gcf,cmap());
clim([min(X),max(X)])
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);
colorbar('EastOutside');

end

function [handle1, handle2, handle3] = drawScatterYraw(Z, X, titlelabel, cmap)
ubound = ceil(max(Z));
lbound = floor(min(Z));
geqone = X >= 1;
leqnone = X <= -1;
other = (X >= -1) & (X <= 1);
colormap(gcf,cmap());
handle2 = scatter(Z(other,1), Z(other,2), 10, X(other), 'filled', 'o');
hold on
handle1 = scatter(Z(geqone,1), Z(geqone,2), 20, X(geqone), '^');
handle3 = scatter(Z(leqnone,1), Z(leqnone,2), 20, X(leqnone), 'v');
hold off
clim([min(X),max(X)])
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);
colorbar('EastOutside');

end
% =========================================================================
function handle = drawScatter(Z, X, titlelabel, cmap)

ubound = ceil(max(Z));
lbound = floor(min(Z));
handle = scatter(Z(:,1), Z(:,2), 8, X, 'filled');
colormap(gcf,cmap());
clim([0,1])
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);
colorbar('EastOutside');

end
% =========================================================================
function drawPortfolioSelections(Z, P, algolabels, titlelabel, cmap)

ubound = ceil(max(Z));
lbound = floor(min(Z));
nalgos = length(algolabels);
algolbls = cell(1,nalgos+1);
h = zeros(1,nalgos+1);
isworthy = sum(bsxfun(@eq, P, 0:nalgos))~=0;
clr = cmap(nalgos+1);
clr = clr([2,1,3:nalgos+1],:);
markers = ['o',repmat('o',1,nalgos)];
clf
hold on
%clr = flipud(lines(nalgos+1));
for i=0:nalgos
    if isworthy(i+1)
        scatter(Z(P==i,1), Z(P==i,2), 8, clr(i+1,:), 'filled', markers(i+1))
        % line(Z(P==i,1), Z(P==i,2), 'LineStyle', 'none', ...
        %                            'Marker', markers(i+1), ...
        %                            'Color', clr(i+1,:), ...
        %                            'MarkerFaceColor', clr(i+1,:), ...
        %                            'MarkerSize', 4);
    end
    h(i+1) = patch([0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
    if i==0
        algolbls{i+1} = 'None';
    else
        algolbls{i+1} = strrep(algolabels{i},'_',' ');
    end
end
hold off
%colormap(gcf,cmap());
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
%legend(h(isworthy), algolbls(isworthy), 'Location', 'NorthEastOutside');
legend(h, algolbls, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================
function h = drawPortfolioFootprint(Z, best, P, algolabels)

% Color definitions
ubound = ceil(max(Z));
lbound = floor(min(Z));
nalgos = length(algolabels);
algolbls = cell(1,nalgos+1);
isworthy = sum(bsxfun(@eq, P, 0:nalgos))~=0;
clr = flipud(lines(nalgos+1));
h = zeros(1,nalgos+1);
for i=0:nalgos
    if ~isworthy(i+1)
        continue;
    end
    line(Z(P==i,1), Z(P==i,2), 'LineStyle', 'none', ...
                               'Marker', '.', ...
                               'Color', clr(i+1,:), ...
                               'MarkerFaceColor', clr(i+1,:), ...
                               'MarkerSize', 4);
    h(i+1) = patch([0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
    if i==0
        algolbls{i+1} = 'None';
    else
        drawFootprint(best{i}, clr(i+1,:), 0.3);
        algolbls{i+1} = strrep(algolabels{i},'_',' ');
    end
end
xlabel('z_{1}'); ylabel('z_{2}'); title('Portfolio footprints');
legend(h(isworthy), algolbls(isworthy), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================
function h = drawGoodBadFootprint(Z, good, Ybin, titlelabel)

ubound = ceil(max(Z));
lbound = floor(min(Z));
orange = [1.0 0.6471 0.0];
blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if any(~Ybin)
    % drawFootprint(bad, orange, 0.2);
    line(Z(~Ybin,1), Z(~Ybin,2), 'LineStyle', 'none', ...
                                 'Marker', '.', ...
                                 'Color', orange, ...
                                 'MarkerFaceColor', orange, ...
                                 'MarkerSize', 4);
    h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
end
if any(Ybin)
    line(Z(Ybin,1), Z(Ybin,2), 'LineStyle', 'none', ...
                               'Marker', '.', ...
                               'Color', blue, ...
                               'MarkerFaceColor', blue, ...
                               'MarkerSize', 4);
    h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
    drawFootprint(good, blue, 0.3);
end
xlabel('z_{1}'); ylabel('z_{2}'); title([titlelabel ' Footprints']);
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================
function handle = drawFootprint(footprint, color, alpha)
% 
hold on;
if isempty(footprint) || isempty(footprint.polygon)
    handle = patch([0 0],[0 0], color, 'EdgeColor','none');
    return
end

handle = plot(footprint.polygon,'FaceColor', color, 'EdgeColor','none', 'FaceAlpha', alpha);
hold off;

end
% =========================================================================
function h = drawBinaryPerformance(Z, Ybin, titlelabel, cmap)

ubound = ceil(max(Z));
lbound = floor(min(Z));
%orange = [1.0 0.6471 0.0];
%blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);

clrs = cmap(2);
orange = clrs(2,:);
blue = clrs(1,:);

if any(~Ybin)
    h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
    line(Z(~Ybin,1), Z(~Ybin,2), 'LineStyle', 'none', ...
                                 'Marker', 'o', ...
                                 'Color', orange, ...
                                 'MarkerFaceColor', orange, ...
                                 'MarkerSize', 2);
end
if any(Ybin)
    h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
    line(Z(Ybin,1), Z(Ybin,2), 'LineStyle', 'none', ...
                               'Marker', 'o', ...
                               'Color', blue, ...
                               'MarkerFaceColor', blue, ...
                               'MarkerSize', 2);
end
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================