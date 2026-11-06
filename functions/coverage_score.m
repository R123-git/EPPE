function [F,CN,CW,U]=coverage_score(S,p)
indices=unique([p.existing(:);S(:)]);
U=any(p.coverage(indices,:),1);
CN=mean(U); CW=sum(double(p.demandWeights(U))); F=.5*(CN+CW);
end
