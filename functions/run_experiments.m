function records=run_experiments(methods,runs,c)
validateattributes(runs,{'numeric'},{'scalar','integer','>=',1,'<=',20});
out=fullfile(c.root,'results','runs');if ~isfolder(out),mkdir(out);end
records={};
for ni=1:2
 for ci=1:2
  net=c.networks{ni};condition=c.conditions{ni}{ci};p=load_instance(c,net,condition);
  for run=1:runs
   seed=c.seedBase+ni*100000+ci*1000+run;
   for j=1:numel(methods)
    method=methods{j};file=fullfile(out,sprintf('%s_%s_%s_%02d.mat',net,condition,method,run));
    if isfile(file)
     a=load(file,'r');r=a.r;assert(r.seed==seed&&strcmp(r.instanceHash,p.signature));
    else
     r=run_method(method,p,c,seed);r.network=net;r.condition=condition;r.run=run;r.instanceHash=p.signature;
     save(file,'r');
    end
    v=layout_metrics(r.S,p);
    assert(r.queries==c.FE.(net)&&r.actualScores+r.cacheHits==r.queries);
    assert(r.feasible&&max(abs([r.F r.Q]-[v.F v.Q]))<1e-10);
    records{end+1}=r; %#ok<AGROW>
    fprintf('%s %s %-14s run %02d: Q=%.6f%% F=%.6f%%\n',net,condition,method,run,100*r.Q,100*r.F);
   end
  end
 end
end
end
