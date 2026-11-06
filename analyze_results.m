function [summary,tests]=analyze_results(source)
% Summarize archived or newly executed runs; incomplete families are not tested.
if nargin<1,source='archive';end;c=config();records=load_records(c,source);
rows={};
for j=1:numel(records)
 r=records{j};rows(end+1,:)={string(r.network),string(r.condition),string(r.method),r.run,100*r.Q,100*r.F,r.feasible,r.seconds,r.queries,r.actualScores,r.cacheHits,r.proxyCalls}; %#ok<AGROW>
end
raw=cell2table(rows,'VariableNames',{'Network','Condition','Method','Run','Q','F','Feasible','ArchivedSeconds','FE','ActualScores','CacheHits','ProxyCalls'});
summary={};tests={};rng(c.statsSeed,'twister');ix=randi(20,20,7000);
for ni=1:2,for ci=1:2
 net=string(c.networks{ni});cond=string(c.conditions{ni}{ci});t=raw(raw.Network==net&raw.Condition==cond,:);
 ep=sortrows(t(t.Method=="EPPE",:),'Run');
 for method=unique(t.Method,'stable')'
  b=sortrows(t(t.Method==method,:),'Run');
  summary(end+1,:)={net,cond,method,height(b),mean(b.Q),std(b.Q),min(b.Q),max(b.Q),mean(b.F),100*mean(b.Feasible),mean(b.ArchivedSeconds),mean(b.FE),mean(b.ProxyCalls)}; %#ok<AGROW>
  if height(ep)~=20||height(b)~=20||method=="EPPE",continue;end
  assert(isequal(ep.Run,b.Run));delta=ep.Q-b.Q;mu=mean(delta);draws=mean(delta(ix),1);ci95=prctile(draws(1:2000),[2.5 97.5]);pv=(1+sum(abs(draws(2001:end)-mu)>=abs(mu)-1e-12))/5001;
  if ismember(method,string(c.main)),add('Main');end
  if ismember(method,["PPE","PPE_M1","PPE_M2","RandomExchange"]),add('Ablation');end
 end
end,end
summary=cell2table(summary,'VariableNames',{'Network','Condition','Method','Runs','MeanQ','SDQ','WorstQ','BestQ','MeanF','FeasiblePercent','ArchivedSeconds','MeanFE','MeanProxyCalls'});
if isempty(tests)
 tests=table();
else
 tests=cell2table(tests,'VariableNames',{'Network','Condition','Family','Reference','DeltaQpp','CILow','CIHigh','P','HolmP','ReferenceFeasiblePercent'});
 for family=["Main","Ablation"]
  use=tests.Family==family;expected=20;if family=="Ablation",expected=16;end
  if sum(use)==expected
   pv=tests.P(use);[v,o]=sort(pv);m=numel(v);adj=zeros(size(v));adj(o)=min(1,cummax(v.*(m:-1:1)'));tests.HolmP(use)=adj;
  end
 end
end
out=fullfile(c.root,'results',source);if ~isfolder(out),mkdir(out);end
writetable(raw,fullfile(out,'raw.csv'));writetable(summary,fullfile(out,'summary.csv'));
if ~isempty(tests),writetable(tests,fullfile(out,'tests.csv'));end
disp(summary);
 function add(family)
  tests(end+1,:)={net,cond,string(family),method,mu,ci95(1),ci95(2),pv,NaN,100*mean(b.Feasible)};
 end
end
