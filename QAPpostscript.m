suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
supplabels = supp.Properties.VariableNames;
issubsource = strcmpi(supplabels,'subsource');
subS = categorical(supp{:,issubsource});
%figure
clf
%set(gcf, 'Position',  [0, 100, 800, 800])
drawSubSources(model.pilot.Z, subS, model.data.S, rootdir);
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
%set(gcf, 'Position',  [0, 100, 800, 800])
Z = model.pilot.Z;
BMAbad = (model.data.Ybin(:,2) == 0);
MMASbad = (model.data.Ybin(:,3) == 0);
%handle = scatter(Z(:,1), Z(:,2), 8, X, 'filled');
hold on
scatter(Z(:,1), Z(:,2),10,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
scatter(Z(BMAbad,1), Z(BMAbad,2),120,'Marker','x','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
scatter(Z(MMASbad,1), Z(MMASbad,2),100,'Marker','+','MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1]);
hold off
%caxis([0,1])
xlabel('z_{1}'); ylabel('z_{2}');
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([-4 4 -4 4]);
print(gcf,'-depsc',[rootdir 'perfBMA_MMAS.eps']);
%colorbar('EastOutside');

function handle = drawSubSources(Z, subS, supS, rootdir)

clf
ubound = ceil(max(Z));
lbound = floor(min(Z));
sourcelabels = cellstr(unique(subS));
nsources = length(sourcelabels);
clrs = flipud(colorcube(nsources+4));
clrs = clrs(2:end-3,:);
clrs = clrs([(0:4)*4+1,(0:4)*4+2,(0:4)*4+3,(0:4)*4+4],:);
handle = zeros(nsources,1);

supercats = categories(supS);
%symbols = {'o','x','s', '^', '+', "pentagram", '<'};
symbols = {'+','s','*', '^', 'o', "pentagram", 'hexagram'};
sizes = [5,5,5,5,5,5,5];

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