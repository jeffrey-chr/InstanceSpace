rootdir = '.\QAPdata\';
model = load('.\QAPdata\model.mat');
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);

easy = [];
hard = [];

for i = 1:length(supp.subsource)
    if strcmp(supp.subsource{i}, 'qapsat-gen-easy') == 1
        easy = [easy,i];
    elseif strcmp(supp.subsource{i}, 'qapsat-gen-hard') == 1
        hard = [hard,i];
    end
end

nf = size(model.data.Xraw,2);
allfdata = [];

for f = 1:nf
    fdata = model.data.Xraw(:,f);
    mx = max(max(abs(fdata)),0.0001);
    allfdata = [allfdata, fdata(easy, :)./mx, fdata(hard, :)./mx];
end

boxplot(allfdata, 'PlotStyle', 'compact')