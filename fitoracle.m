function out = fitoracle(Z,Y,opts)

% global params
Y = double(Y)+1;
nalgos = size(Y,2);
out.paramgrid = sortrows(2.^(opts.maxcvgrid.*lhsdesign(opts.cvgrid,2) + ...
                             opts.mincvgrid));  % Cross-validation grid
out.cvmcr = NaN.*ones(opts.cvgrid,nalgos);
out.paramidx = NaN.*ones(1,nalgos);

for i=1:nalgos
    for j=1:opts.cvgrid
        out.cvmcr(j,i) = crossval('mcr', Z, Y(:,i),...
                                  'Kfold', opts.cvfolds,...
                                  'Options', statset('UseParallel',true),...
                                  'Predfun',@(xtrain,ytrain,xtest) svmwrap(xtrain,...
                                                                           ytrain,...
                                                                           xtest,...
                                                                           out.paramgrid(j,:)));
    end
    [~,out.paramidx(i)] = min(out.cvmcr(:,i));
    disp(['    ' num2str(i) ' out of ' num2str(nalgos) ' models have been fitted.']);
end