newdir = './QAPdata/';
newfeatfile = [newdir 'metadata_supp.csv'];
newfeattable = readtable(newfeatfile);
newfeat = table2array(newfeattable(:,2:end-3));
newalg = table2array(newfeattable(:,end-2:end-1));
newfeatraw = newfeat;

fullalg = newalg(:,1) - newalg(:,2);

edgesDre = 8*9+10*7;
edgesHyp = 4*3*3*3*2;

ntriDre = zeros(30,1);
ntriHyp = zeros(30,1);
ntriDre(1) = edgesDre;
ntriHyp(1) = edgesHyp;
for i = 2:30
    tmp = (i-1)/29;
    ntriDre(i) = floor(edgesDre + tmp*(edgesDre*3));
    ntriHyp(i) = floor(edgesHyp + tmp*(edgesHyp*3));
end

nSeeds = 6:35;
nSquares = ceil(5.5:0.5:20);

figure(1)
scatter(ntriDre, fullalg(61:90))
hold on
fplot(@(x) 0, [ntriDre(1) ntriDre(end)], 'k')
fplot(@(x) 1, [ntriDre(1) ntriDre(end)], 'b')
fplot(@(x) -1, [ntriDre(1) ntriDre(end)], 'r')
hold off
ylim([-2 2]);
title("N cycles vs performance, Drexx distances PROBABLY WRONG")

figure(2)
scatter(ntriHyp, fullalg(91:120))
hold on
fplot(@(x) 0, [ntriHyp(1) ntriHyp(end)], 'k')
fplot(@(x) 1, [ntriHyp(1) ntriHyp(end)], 'b')
fplot(@(x) -1, [ntriHyp(1) ntriHyp(end)], 'r')
hold off
ylim([-2 2]);
title("N cycles vs performance, Hcube distances PROBABLY WRONG")

if false
figure(5)
sp = newfeattable.feature_FlowSparsity(1:30);
scatter(sp, fullalg(61:90))
hold on
fplot(@(x) 0, [min(sp) max(sp)], 'k')
fplot(@(x) 1, [min(sp) max(sp)], 'b')
fplot(@(x) -1, [min(sp) max(sp)], 'r')
hold off
ylim([-2 2]);
title("Sparsity of tcycle instance vs performance, Drexx distances")

figure(6)
sp = newfeattable.feature_FlowSparsity(31:60);
scatter(sp, fullalg(91:120))
hold on
fplot(@(x) 0, [min(sp) max(sp)], 'k')
fplot(@(x) 1, [min(sp) max(sp)], 'b')
fplot(@(x) -1, [min(sp) max(sp)], 'r')
hold off
ylim([-2 2]);
title("Sparsity of tcycle instance vs performance, Hcube distances")
end

figure(3)
scatter(nSeeds, fullalg(121:150))
hold on
fplot(@(x) 0, [nSeeds(1) nSeeds(end)], 'k')
fplot(@(x) 1, [nSeeds(1) nSeeds(end)], 'b')
fplot(@(x) -1, [nSeeds(1) nSeeds(end)], 'r')
hold off
ylim([-2 2]);
title("N webseeds vs performance, Drexx distances PROBABLY WRONG")

figure(4)
scatter(nSeeds, fullalg(151:180))
hold on
fplot(@(x) 0, [nSeeds(1) nSeeds(end)], 'k')
fplot(@(x) 1, [nSeeds(1) nSeeds(end)], 'b')
fplot(@(x) -1, [nSeeds(1) nSeeds(end)], 'r')
hold off
ylim([-2 2]);
title("N webseeds vs performance, Hcube distances PROBABLY WRONG")

figure(5)
scatter(nSquares, fullalg(1:30))
hold on
fplot(@(x) 0, [nSquares(1) nSquares(end)], 'k')
fplot(@(x) 1, [nSquares(1) nSquares(end)], 'b')
fplot(@(x) -1, [nSquares(1) nSquares(end)], 'r')
hold off
ylim([-2 2]);
title("N squares vs performance, Drexx distances PROBABLY WRONG")

figure(6)
scatter(nSquares, fullalg(31:60))
hold on
fplot(@(x) 0, [nSquares(1) nSquares(end)], 'k')
fplot(@(x) 1, [nSquares(1) nSquares(end)], 'b')
fplot(@(x) -1, [nSquares(1) nSquares(end)], 'r')
hold off
ylim([-2 2]);
title("N squares vs performance, Hcube distances PROBABLY WRONG")