function [Q,F,loss,worst]=robust_score(S,p)
T=sort([p.existing S]); C=double(p.coverage(T,:)); count=sum(C,1); w=.5/numel(p.ids)+.5*p.demandWeights;
F=sum(w(count>0)); losses=C*(w.*(count==1))'; [loss,worst]=max(losses); Q=F-loss;
end
