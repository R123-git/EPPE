function report=diagnostics(mode,runs,seconds)
% Optional analyses used in the article; main searches remain separate.
if nargin<1,mode='verify';end;if nargin<2,runs=20;end;if nargin<3,seconds=60;end
c=config();out=fullfile(c.root,'results','diagnostics');if ~isfolder(out),mkdir(out);end
validateattributes(runs,{'numeric'},{'scalar','integer','>=',1,'<=',20});
rows={};
switch mode
 case {'verify','bounds'}
  for ni=1:2,for ci=1:2
   net=c.networks{ni};cond=c.conditions{ni}{ci};p=load_instance(c,net,cond);
   z=load(fullfile(c.root,'data','diagnostics',[net '_' cond '_final.mat']),'results');
   for k=1:2
    r=z.results{k};
    if strcmp(mode,'bounds'),r=robust_milp(p,r.floor,r.S,seconds);save(fullfile(out,sprintf('%s_%s_L%d.mat',net,cond,2*(k-1))),'r');end
    v=layout_metrics(r.S,p);
    assert(max(abs([r.F r.lowerQ]-[v.F v.Q]))<1e-8&&v.F>=r.floor-1e-9&&r.lowerQ<=r.upperQ+1e-8);
    rows(end+1,:)={string(net),string(cond),2*(k-1),100*v.F,100*v.Q,100*r.upperQ,r.certified}; %#ok<AGROW>
   end
  end,end
  report=cell2table(rows,'VariableNames',{'Network','Condition','LossAllowancePP','F','LowerQ','UpperQ','Certified'});
 case {'budget','threshold'}
  files=dir(fullfile(c.root,'data','diagnostics','Modena_*_B*.mat'));methods={'EPPE','PPE'};baseSeed=2720110101;
  if strcmp(mode,'threshold'),files=dir(fullfile(c.root,'data','diagnostics','Ctown_*_T*.mat'));methods={'EPPE','RandomExchange'};baseSeed=2730110101;end
  for pi=1:numel(files)
   a=load(fullfile(files(pi).folder,files(pi).name),'p');p=a.p;p.objective='QFloor';
   for run=1:runs,for m=1:numel(methods)
    r=run_method(methods{m},p,c,baseSeed+pi*1000+run);r.run=run;r.problem=erase(files(pi).name,'.mat');r.network=p.network;r.condition=p.condition;
    save(fullfile(out,sprintf('%s_%s_%02d.mat',r.problem,r.method,run)),'r');
    rows(end+1,:)={string(r.problem),string(r.method),run,100*r.F,100*r.Q,r.queries,r.actualScores,r.cacheHits,r.proxyCalls}; %#ok<AGROW>
   end,end
  end
  report=cell2table(rows,'VariableNames',{'Problem','Method','Run','F','Q','FE','ActualScores','CacheHits','ProxyCalls'});
 case 'sensitivity'
  a=load(fullfile(c.root,'data','archive','records.mat'),'records');
  for ni=1:2,for ci=1:2
   net=c.networks{ni};cond=c.conditions{ni}{ci};nominal=load_instance(c,net,cond);
   selected=a.records(cellfun(@(r)strcmp(r.network,net)&&strcmp(r.condition,cond)&&ismember(r.method,c.main)&&r.run<=5,a.records));
   for panel={'Threshold','Demand'}
    values=[.75 .85];if strcmp(panel{1},'Demand'),values=[.95 1.05];end
    for value=values
     p=nominal;
     if strcmp(panel{1},'Threshold')
      p.coverage=p.distance<=2000&p.H<=20&p.R>=value&(p.validResponse&p.validResponse');
     else
      z=load(fullfile(c.root,'data','diagnostics',sprintf('%s_scale%03d.mat',net,round(100*value))),'p');p.coverage=z.p.coverage;p.demandWeights=z.p.demandWeights;
     end
     for j=1:numel(selected)
      r=selected{j};v=layout_metrics(r.S,p);
      rows(end+1,:)={string(net),string(cond),string(r.method),r.run,string(panel{1}),value,100*v.F,100*v.Q,v.F>=nominal.floor-1e-10}; %#ok<AGROW>
     end
    end
   end
  end,end
  report=cell2table(rows,'VariableNames',{'Network','Condition','Method','Run','Panel','Value','F','Q','PassNominalFloor'});
 otherwise,error('Use verify, bounds, budget, threshold or sensitivity.');
end
writetable(report,fullfile(out,[mode '.csv']));
fprintf('%s: saved %d records under results/diagnostics.\n',mode,height(report));
if height(report)<=8,disp(report);end
end
