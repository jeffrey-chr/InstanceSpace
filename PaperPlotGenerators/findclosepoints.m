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

target = [-1, -1];

%target = [1, -2];
%type = "Hypercube";
%type = "Terminal";

% target = [-2, -1];
% type = "Palubeckis";
% 
%target = [-2, -2];
%target = [2, 2];
%type = "RealLifeLike";
%type = "RealLife"

%target = [-1.5, -0.7];
%target = [1,0];
%type = "Other";

%target = [-0.32, -0.75];
%type = "RealLifeLike";

%target = [3,0];
%type = "FlowCluster";

%target = [-2,-2];
%type = "RandomUniform";

%target = [-2, -1];
%type = "ManhattanDist";

% target = [-2, -1];
% type = "Other";


%lipa40b, weird on left side

%sko49, manhattan on left side
%dre56, drezner on right
%stf60es2, real-like on right
%stf80er3 or stf60er1, real-like on left
%hyp64_3, hypercube in bottom right
%term45_9, terminal in bottom right

%cand = find(model.data.S == type);
cand = 1:length(model.data.S);

names = model.data.instlabels(cand);
points = model.pilot.Z(cand,:);
dists = zeros(size(points,1),1);
for i = 1:length(points)
    dists(i) = norm(points(i,:) - target);
end

[~,II] = sort(dists,'ascend');

for i = 1:30
    j = II(i);
    fprintf("%s: %f (%f, %f)\n", names{j}, dists(j), points(j,1), points(j,2));
end