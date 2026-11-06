function [p,certificate]=normal_solution(p,B)
p.B=B; p.k=B-numel(p.existing); assert(p.k>=1&&p.k<=numel(p.available)); n=numel(p.available); m=numel(p.ids); w=.5/m+.5*p.demandWeights;
existing=double(any(p.coverage(p.existing,:),1))'; C=sparse(double(p.coverage(p.available,:))');
options=optimoptions('intlinprog','Display','off','MaxTime',60,'RelativeGapTolerance',1e-9,'AbsoluteGapTolerance',1e-9,'ConstraintTolerance',1e-9);
timer=tic; [x,value,flag,out]=intlinprog([zeros(n,1);-w(:)],1:n,[-C speye(m)],existing,[ones(1,n) zeros(1,m)],p.k,zeros(n+m,1),ones(n+m,1),options); elapsed=toc(timer);
assert(flag==1&&out.absolutegap<=1e-8,'Normal F needs a certificate'); S=sort(p.available(x(1:n)>.5)); [~,F]=robust_score(S,p); assert(abs(F+value)<1e-8);
p.Fstar=F; p.floor=F-.02; p.floorTolerance=1e-10; p.initialSet=S;
[~,pos]=ismember(S,p.available); rest=setdiff(1:n,pos,'stable'); p.initialKeys=zeros(1,n); p.initialKeys([pos rest])=((1:n)-.5)/n;
p.F0=coverage_score([],p); p.signature=study_hash(struct('base',p.signature,'budget',B,'coverage',p.coverage,'Fstar',F,'S',S));
certificate=struct('flag',flag,'out',out,'seconds',elapsed,'Fstar',F,'S',S);
end
