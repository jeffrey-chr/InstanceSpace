rootdir = '.\QAPdata\';
model = load('.\QAPdata\model.mat');
suppfile = [rootdir 'suppdata.csv'];
supp = readtable(suppfile);
unftfile = [rootdir 'UNSCALED_featdata.csv'];
unft = readtable(unftfile);

% match instance labels

[names, mnI] = sort(model.data.instlabels);
[~, ufI] = sort(unft.instances);
[~, spI] = sort(supp.instances);

if ~all(strcmp(model.data.instlabels(mnI), supp.instances(spI)))
    error("Instance names do not match");
end

% 

nI = size(model.pilot.Z,1);
guesslen = 2000;

%distances = Inf*ones(nI);
%distances_short = sparse(nI, nI);
distances_i = -ones(guesslen,1);
distances_j = -ones(guesslen,1);
distances_d = -ones(guesslen,1);
count = 0;

gridsize = 0.25;
rounded = floor(model.pilot.Z / gridsize) * gridsize;
minis = min(rounded);
maxes = max(rounded);
griddim = (maxes - minis) / gridsize + 1;
gridid = (rounded - minis) / gridsize + 1;

alloc = cell(griddim);
for i = 1:nI
    alloc{gridid(i,1),gridid(i,2)} = [ alloc{gridid(i,1),gridid(i,2)}; i];
end



tic
%dirs = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1];
dirs = [-1,-1; -1,0; -1,1; 0,-1]; % only check forwards
for g1 = 1:griddim(1)
    for g2 = 1:griddim(2)

        % check same grid box
        for ig = 1:size(alloc{g1,g2},1)
            for ia = (ig+1):size(alloc{g1,g2},1)
                i = alloc{g1,g2}(ig);
                j = alloc{g1,g2}(ia);
                count = count + 1;
                if count > length(distances_d)
                    distances_i = [distances_i; -ones(count-1,1)]; %#ok<AGROW> 
                    distances_j = [distances_j; -ones(count-1,1)]; %#ok<AGROW> 
                    distances_d = [distances_d; -ones(count-1,1)]; %#ok<AGROW> 
                end
                distances_i(count) = i;
                distances_j(count) = j;
                distances_d(count) = norm(model.pilot.Z(i,:) - model.pilot.Z(j,:));
            end
        end
        
        % check adjacent boxes
        for d = 1:size(dirs,1)
            adj1 = g1 + dirs(d,1);
            adj2 = g2 + dirs(d,2);
            if adj1 > 0.5 && adj1 < griddim(1) + 0.5 && adj2 > 0.5 && adj2 < griddim(2) + 0.5
                for ig = 1:size(alloc{g1,g2},1)
                    for ia = 1:size(alloc{adj1,adj2},1)
                        i = alloc{g1,g2}(ig);
                        j = alloc{adj1,adj2}(ia);
                        count = count + 1;
                        if count > length(distances_d)
                            distances_i = [distances_i; -ones(count-1,1)]; %#ok<AGROW> 
                            distances_j = [distances_j; -ones(count-1,1)]; %#ok<AGROW> 
                            distances_d = [distances_d; -ones(count-1,1)]; %#ok<AGROW> 
                        end
                        distances_i(count) = i;
                        distances_j(count) = j;
                        distances_d(count) = norm(model.pilot.Z(i,:) - model.pilot.Z(j,:));
                    end
                end
            end
        end
    end
end
distances_i = distances_i(1:count);
distances_j = distances_j(1:count);
distances_d = distances_d(1:count);
toc

[distances_d, dorder] = sort(distances_d);
distances_i = distances_i(dorder);
distances_j = distances_j(dorder);

algo = 1;
for pid = 1:length(distances_d)
    inst1 = distances_i(pid);
    inst2 = distances_j(pid);
    good1 = model.data.Ybin(inst1,algo);
    good2 = model.data.Ybin(inst2,algo);
    if good1 ~= good2
        fprintf("%s [%0.2f,%0.2f] (%0.4f) vs. %s [%0.2f,%0.2f] (%0.4f)\n",model.data.instlabels{inst1}, model.pilot.Z(inst1,1), model.pilot.Z(inst1,2), model.data.Yraw(inst1,algo), model.data.instlabels{inst2}, model.pilot.Z(inst2,1), model.pilot.Z(inst2,2), model.data.Yraw(inst2,algo))
    end
end


% tic
% for i = 1:nI
%     for j = 1:nI
%         distances(i,j) = norm(model.pilot.Z(i,:) - model.pilot.Z(j,:));
%     end
% end
% toc

% for i = 1:nI
%     for j = 1:nIite(distances2(i,j))
%             i
%             j
%             error();
%         end
%     end
% end
%         if distances(i,j) ~= distances2(i,j) && isfin