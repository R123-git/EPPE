function r=robust_milp(p,floorValue,warmSet,limit)
n=numel(p.available); m=numel(p.ids); B=p.B; w=.5/m+.5*p.demandWeights; C=sparse(double(p.coverage(p.available,:))'); a=double(sum(p.coverage(p.existing,:),1))';
cols=n+2*m+1; iy=n+(1:m); iz=n+m+(1:m); iq=cols; sensors=sort([p.existing p.available]); ns=numel(sensors);
A=spalloc(4*m+2+ns,cols,4*nnz(C)+8*m+ns*(2*m+2)); rhs=zeros(size(A,1),1); I=speye(m);
A(1:m,1:n)=-C; A(1:m,iy)=I; rhs(1:m)=a;
A(m+(1:m),1:n)=C; A(m+(1:m),iy)=-B*I; rhs(m+(1:m))=-a;
A(2*m+(1:m),1:n)=-C; A(2*m+(1:m),iz)=2*I; rhs(2*m+(1:m))=a;
A(3*m+(1:m),1:n)=C; A(3*m+(1:m),iz)=-(B-1)*I; rhs(3*m+(1:m))=1-a;
A(4*m+1,iy)=-w; rhs(4*m+1)=-floorValue; A(4*m+2,iy)=-w; A(4*m+2,iq)=1;
for j=1:ns
    row=4*m+2+j; covered=logical(p.coverage(sensors(j),:)); A(row,iy(~covered))=-w(~covered); A(row,iz(covered))=-w(covered); A(row,iq)=1;
    pos=find(p.available==sensors(j)); if ~isempty(pos), bigM=sum(w(covered)); A(row,pos)=bigM; rhs(row)=bigM; end
end
lower=zeros(cols,1); upper=ones(cols,1);
if floorValue>sum(w)-min(w)+1e-10, lower(iy)=1; end
T=sort([p.existing warmSet]); count=sum(p.coverage(T,:),1); x0=zeros(cols,1); x0(ismember(p.available,warmSet))=1; x0(iy)=count>0; x0(iz)=count>1; x0(iq)=robust_score(warmSet,p);
assert(sum(w(count>0))>=floorValue-1e-10&&numel(warmSet)==p.k&&all(A*x0<=rhs+1e-8));
eq=sparse(1,1:n,1,1,cols); objective=zeros(cols,1); objective(iq)=-1;
options=optimoptions('intlinprog','Display','off','MaxTime',limit,'RelativeGapTolerance',1e-9,'AbsoluteGapTolerance',1e-9,'ConstraintTolerance',1e-9);
timer=tic; [x,value,flag,out]=intlinprog(objective,1:n+2*m,A,rhs,eq,p.k,lower,upper,x0,options); elapsed=toc(timer);
fallback=isempty(x); if fallback, S=warmSet; upperQ=1; auxiliaryQ=NaN; else, S=sort(p.available(x(1:n)>.5)); auxiliaryQ=-value; upperQ=min(1,-value+out.absolutegap); end
[v,failures,failedSensors]=layout_metrics(S,p); [Q,F]=robust_score(S,p);
assert(numel(unique(S))==p.k&&all(ismember(S,p.available))&&F>=floorValue-1e-9&&max(abs([Q F]-[v.Q v.F]))<1e-10);
assert(upperQ>=Q-1e-7,'Bound must contain independently checked feasible Q'); upperQ=max(Q,upperQ);
r=struct('S',S,'F',F,'Q',Q,'lowerQ',Q,'upperQ',upperQ,'flag',flag,'out',out,'seconds',elapsed,'limit',limit,'floor',floorValue,'certified',flag==1&&upperQ-Q<1e-7,'fallback',fallback,'auxiliaryQ',auxiliaryQ,'fullyCovered',all(any(p.coverage([p.existing S],:),1)),'failures',failures,'failedSensors',p.ids(failedSensors),'variables',cols,'integerVariables',n+2*m,'constraints',size(A,1));
end
