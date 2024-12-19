clear;

outputdir = '.\output\';

f = gcf;
f.Position = [50 750 800 600];

cmap = @copper;

newdir = '..\QAP_separateexp\';
fcfeatfile = [newdir 'metadata_fc.csv'];
fcfeattable = readtable(fcfeatfile);
fcfeat = table2array(fcfeattable(:,2:end-3));
fcalg = table2array(fcfeattable(:,end-2:end-1));
fcsuppfile = [newdir 'suppdata_fc.csv'];
fcsupptable = readtable(fcsuppfile);
fcsubS = categorical(fcsupptable.subsource);

fcfullalg = fcalg(:,2) - fcalg(:,1);

nSeeds = 16:45;
nSquares = ceil(5.5:0.5:20);
nTri = 200:25:925;

bk = repmat([0,0,0],30,1);

figure(1)
hc1 = scatter(nTri, fcfullalg(fcsubS == 'flowcluster-dhyper-fcycle'), 25, bk, '*')
hold on
hc2 = scatter(nTri, fcfullalg(fcsubS == 'flowcluster-ddrez-fcycle') , 25, bk, 's')
hc3 = fplot(@(x) 0*x+1, [175 950], '--k')
hc4 = fplot(@(x) 0*x-1, [175 950], ':k')
hc5 = fplot(@(x) 0*x, [175 950], 'k')
hold off
ylim([-2 2]);
xlim([175 950]);
title("Triangle flow instances")
xlabel("Number of triangles");
ylabel("Algorithm performance")
legend([hc1,hc2,hc3,hc4,hc5], ["Hcube distances", "Drexx distances", "BMA stronger", "MMAS stronger", "Algorithms even"], 'Location', 'SouthOutside', 'NumColumns', 3);
print(gcf,'-dpng',[outputdir 'ext_perf_tcycle.png']);
print(gcf,'-depsc',[outputdir 'ext_perf_tcycle.eps']);

hs1 = scatter(nSquares, fcfullalg(fcsubS == 'flowcluster-dhyper-fsquare'), 25, bk, '*')
hold on
hs2 = scatter(nSquares, fcfullalg(fcsubS == 'flowcluster-ddrez-fsquare'), 25, bk, 's')
hs3 = fplot(@(x) 0*x+1, [5 21], '--k')
hs4 = fplot(@(x) 0*x-1, [5 21], ':k')
hs5 = fplot(@(x) 0*x, [5 21], 'k')
hold off
xlim([5 21]);
ylim([-2 2]);
xlabel("Number of squares");
ylabel("Algorithm performance")
legend([hs1,hs2,hs3,hs4,hs5], ["Hcube distances", "Drexx distances", "BMA stronger", "MMAS stronger", "Algorithms even"], 'Location', 'SouthOutside', 'NumColumns', 3);
title("Square-based flow instances")
print(gcf,'-dpng',[outputdir 'ext_perf_square.png']);
print(gcf,'-depsc',[outputdir 'ext_perf_square.eps']);

ht1 = scatter(nSeeds, fcfullalg(fcsubS == 'flowcluster-dhyper-ftree'), 25, bk, '*')
hold on
ht2 = scatter(nSeeds, fcfullalg(fcsubS == 'flowcluster-ddrez-ftree'), 25, bk, 's')
ht3 = fplot(@(x) 0*x+1, [15 46], '--k')
ht4 = fplot(@(x) 0*x-1, [15 46], ':k')
ht5 = fplot(@(x) 0*x, [15 46], 'k')
hold off
xlim([15 46]);
ylim([-2 2]);
title("Tree-based flow instances")
xlabel("Number of tree seeds");
ylabel("Algorithm performance");
legend([ht1,ht2,ht3,ht4,ht5], ["Hcube distances", "Drexx distances", "BMA stronger", "MMAS stronger", "Algorithms even"], 'Location', 'SouthOutside', 'NumColumns', 3);
print(gcf,'-dpng',[outputdir 'ext_perf_tree.png']);
print(gcf,'-depsc',[outputdir 'ext_perf_tree.eps']);

