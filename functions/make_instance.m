function p=make_instance(base,condition,cfg)
clock=tic; p=base; net=p.network; V=p.eligible; B=cfg.B.(net);
networkIndex=find(strcmp({'Hanoi','Net3','Ctown','Anytown','Modena'},net));
maskNumber=0; type=0;
if startsWith(condition,'R'), type=1; maskNumber=str2double(condition(2:end)); end
if startsWith(condition,'B'), type=2; maskNumber=str2double(condition(2:end)); end
stream=RandStream('mt19937ar','Seed',cfg.maskSeedBase+networkIndex*1000+10*type+maskNumber);
S0=[];
if ~strcmp(condition,'C0')
    [~,i]=max(p.nominalDemand(V)); S0=V(i);
    while numel(S0)<cfg.existing.(net)
        free=setdiff(V,S0); [~,j]=max(min(p.distance(free,S0),[],2)); S0(end+1)=free(j); %#ok<AGROW>
    end
end
free=setdiff(V,S0); forbidden=[];
if startsWith(condition,'R')
    forbidden=free(randperm(stream,numel(free),ceil(cfg.forbiddenFraction*numel(free))));
elseif startsWith(condition,'B')
    centre=free(randi(stream,numel(free))); [~,j]=sort(p.distance(centre,free),'ascend');
    forbidden=free(j(1:ceil(cfg.forbiddenFraction*numel(free))));
elseif ~ismember(condition,{'C0','C1'}), error('未知情景：%s',condition); end
A=setdiff(free,forbidden); k=B-numel(S0); assert(k>0 && numel(A)>=k,'布局不可行');
adj=double(adjacency(p.graph)>0); degree=sum(adj,2); invd=zeros(size(degree));
invd(degree>0)=1./sqrt(degree(degree>0)); lap=eye(size(adj))-invd.*adj.*invd';
parts=conncomp(graph(adj)); labels=zeros(numel(degree),1);
K=max(ceil(sqrt(numel(p.ids))),max(parts)); sizes=accumarray(parts(:),1);
allocation=apportion(sizes,K);
old=rng; restore=onCleanup(@()rng(old)); %#ok<NASGU>
rng(cfg.partitionSeedBase+networkIndex*1000,'twister'); offset=0;
for c=1:numel(sizes)
    rows=find(parts==c); groups=allocation(c);
    if groups==1, local=ones(numel(rows),1);
    else
        [vec,val]=eig(lap(rows,rows),'vector'); [~,order]=sort(val);
        Y=vec(:,order(1:groups)); Y=Y./max(vecnorm(Y,2,2),eps);
        local=kmeans(Y,groups,'Replicates',20,'MaxIter',1000,'EmptyAction','singleton');
    end
    minima=accumarray(local,rows,[],@min); [~,order]=sort(minima); canon=zeros(numel(order),1); canon(order)=1:numel(order);
    labels(rows)=canon(local)+offset; offset=offset+groups;
end
labels=labels(p.junctionIndex); entropy=zeros(1,numel(p.ids));
X=abs(p.deltaPressure(p.validResponse,:)); maximum=max(X,[],'all');
if maximum>0
    edges=linspace(0,maximum,21);
    for j=1:size(X,2)
        q=histcounts(X(:,j),edges)/size(X,1); q=q(q>0); entropy(j)=-sum(q.*log(q));
    end
end
m=min(numel(A),max(k,ceil(cfg.priorityFraction*numel(A))));
groups=unique(labels(A)); capacity=arrayfun(@(g)sum(labels(A)==g),groups);
if numel(groups)>m
    value=arrayfun(@(g)max(entropy(A(labels(A)==g))),groups);
    [~,order]=sort(value,'descend'); groups=groups(order(1:m)); capacity=capacity(order(1:m));
end
quota=apportion(capacity(:),m); P0=[];
for g=1:numel(groups)
    available=A(labels(A)==groups(g)); [~,order]=sort(entropy(available),'descend');
    P0=[P0 available(order(1:quota(g)))]; %#ok<AGROW>
end
p.existing=sort(S0); p.forbidden=sort(forbidden); p.available=A(:)'; p.k=k; p.B=B;
p.P0=sort(P0); p.condition=condition; p.partition=labels; p.entropy=entropy;
p.F0=coverage_score([],p); p.prepareSeconds=toc(clock);
p.signature=study_hash(struct('base',base.signature,'S0',S0,'M',forbidden,'B',B,'P0',p.P0, ...
    'instanceBuilder',study_hash(mfilename('fullpath')+".m",true)));
end
function q=apportion(capacity,total)
q=ones(size(capacity)); assert(total>=numel(q)&&total<=sum(capacity));
while sum(q)<total
    left=total-sum(q); weights=capacity; weights(q>=capacity)=0;
    ideal=left*weights/sum(weights); extra=min(floor(ideal),capacity-q); q=q+extra;
    left=total-sum(q); if left==0, break; end
    remainder=ideal-floor(ideal); remainder(q>=capacity)=-Inf;
    [~,order]=sort(remainder,'descend'); eligible=order(isfinite(remainder(order)));
    take=eligible(1:min(left,numel(eligible))); q(take)=q(take)+1;
end
end
