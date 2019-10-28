function out = findSpaceBounds(X,A,opts)

nfeats = size(X,2);
[rho,pval] = corr(X);
rho = rho.*(pval<opts.pval);

Xbnds = [min(X); max(X)];
idx = de2bi(0:2^nfeats-1)+1;
Xedge = zeros(size(idx,1),nfeats);
flags = false(size(idx,1),1);
for i=1:size(idx,1)
   ind = sub2ind(size(Xbnds),idx(i,:),1:nfeats);
   Xedge(i,:) = Xbnds(ind)';
   for j=1:nfeats
       for k=j+1:nfeats
           % Check for valid points give the correlation trend
           if rho(j,k)>opts.cthres && sign(Xedge(i,j))~=sign(Xedge(i,k))
               flags(i) = true;
           elseif rho(j,k)<-opts.cthres && sign(Xedge(i,j))==sign(Xedge(i,k))
               flags(i) = true;
           end
       end
   end
end
Zedge = Xedge*A';
Kedge = convhull(Zedge(:,1),Zedge(:,2));
out.Zedge = Zedge(Kedge,:);

try
    Xecorr = Xedge(~flags,:);
    Zecorr = Xecorr*A';
    Kecorr = convhull(Zecorr(:,1),Zecorr(:,2));
    out.Zecorr = Zecorr(Kecorr,:);
catch
    disp('  -> The correlation threshold is too strict (your features are highly correlated).')
    disp('  -> Please consider increasing it.');
    out.Zecorr = out.Zedge;
end

