function r=robust_milp_full(p,floorValue,warmSet,limit)
n=numel(p.available); m=numel(p.ids); scenarios=sort([p.existing p.available]); ns=numel(scenarios); w=.5/m+.5*p.demandWeights;
cols=n+(ns+1)*m+1; rows=(ns+1)*m+ns+1; A=spalloc(rows,cols,(ns+1)*(nnz(p.coverage(p.available,:))+3*m)+ns*m); rhs=zeros(rows,1); x0=zeros(cols,1); x0(ismember(p.available,warmSet))=1; T=sort([p.existing warmSet]);
for block=0:ns
    rr=block*m+(1:m); yc=n+block*m+(1:m); C=sparse(double(p.coverage(p.available,:))'); E=p.existing;
    if block>0, failed=scenarios(block); C(:,p.available==failed)=0; E=setdiff(E,failed); U=any(p.coverage(setdiff(T,failed),:),1); else, U=any(p.coverage(T,:),1); end
    A(rr,1:n)=-C; A(rr,yc)=speye(m); rhs(rr)=double(any(p.coverage(E,:),1))'; x0(yc)=double(U);
    row=(ns+1)*m+block+1; A(row,yc)=-w; if block==0, rhs(row)=-floorValue; else, A(row,cols)=1; end
end
x0(cols)=robust_score(warmSet,p); objective=zeros(cols,1); objective(cols)=-1; eq=sparse(1,1:n,1,1,cols); options=optimoptions('intlinprog','Display','off','MaxTime',limit,'RelativeGapTolerance',1e-9,'AbsoluteGapTolerance',1e-9,'ConstraintTolerance',1e-9);
[x,value,flag,out]=intlinprog(objective,1:n,A,rhs,eq,p.k,zeros(cols,1),ones(cols,1),x0,options); S=p.available(x(1:n)>.5); [Q,F]=robust_score(S,p); r=struct('Q',Q,'F',F,'flag',flag,'out',out,'auxQ',-value);
end
