function [v,failures,T]=layout_metrics(S,p)
T=sort([p.existing S]); assert(numel(unique(T))==numel(T)); w=.5/numel(p.ids)+.5*p.demandWeights;
U=any(p.coverage(T,:),1); F=sum(w(U)); CN=mean(U); CW=sum(p.demandWeights(U)); failures=zeros(numel(T),3);
for j=1:numel(T)
    U=any(p.coverage(T([1:j-1 j+1:end]),:),1); failures(j,:)=[sum(w(U)),mean(U),sum(p.demandWeights(U))];
end
[Q,j]=min(failures(:,1)); v=struct('Q',Q,'F',F,'CN',CN,'CW',CW,'QCN',failures(j,2),'QCW',failures(j,3),'worstLoss',F-Q,'worstSensor',p.ids(T(j)));
end
