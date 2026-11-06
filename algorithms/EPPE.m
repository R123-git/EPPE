function result=EPPE(p,c,seed,mode)
% PPE search with optional coverage-guided set exchanges.
if nargin<4, mode='EPPE'; end
useAdd=ismember(mode,{'EPPE','PPE_M1'});
useRemove=ismember(mode,{'EPPE','PPE_M2'});
enableExchange=~strcmp(mode,'PPE');
rng(seed,'twister'); stream=RandStream('mt19937ar','Seed',seed+1000000);
e=RobustEvaluator(p,c.FE.(p.network)); N=c.population; D=numel(p.available);
X=rand(N,D); X(1,:)=p.initialKeys; velocity=zeros(N,D); loss=Inf(N,1); Px=ones(N,1)/N;
Pa=1.1; K=floor(log(N))+1; sigma=.1; TT=sigma; maxGeneration=ceil(e.cap/N);
for i=1:N, loss(i)=e.evaluate(X(i,:)); end
[~,order]=sort(loss); elite=X(order(1:K),:); eliteLoss=loss(order(1:K));
queue=[]; stagnant=0; lastTrigger=-c.stagnation; batch=1; exchangeCount=0;
while e.remaining()
    before=e.bestF; trial=X+velocity; outside=trial<0|trial>1; trial(outside)=rand(sum(outside,'all'),1);
    movement=velocity; childLoss=Inf(N,1); evaluated=false(N,1);
    for i=1:N
        if ~e.remaining(), break; end
        if ismember(i,queue)
            trial(i,:)=exchange(X(i,:),p,c,stream,e,useAdd,useRemove);
            movement(i,:)=trial(i,:)-X(i,:); exchangeCount=exchangeCount+1;
        end
        childLoss(i)=e.evaluate(trial(i,:)); evaluated(i)=true;
    end
    [~,order]=sort(childLoss);
    for k=1:K
        if childLoss(order(k))<=eliteLoss(k), elite(k,:)=trial(order(k),:); eliteLoss(k)=childLoss(order(k)); end
    end
    for i=find(evaluated)'
        distance=vecnorm(elite-X(i,:),2,2); %#ok<NASGU>
        if childLoss(i)<=loss(i)
            X(i,:)=trial(i,:); loss(i)=childLoss(i); Px(i)=Pa*Px(i)*(1-Px(i));
            [~,slot]=min(vecnorm(elite-X(i,:),2,2)); direction=elite(slot,:)-X(i,:);
            A3=zeros(1,D); count=floor(rand*D); if count>0, take=randperm(D,count); A3(take)=2*sigma*randn(1,count); end
            velocity(i,:)=(1-Px(i))*.2*direction+Px(i)*(movement(i,:)+A3);
        else
            if rand<=Px(i), X(i,:)=trial(i,:); loss(i)=childLoss(i); Px(i)=Pa*Px(i)*(1-Px(i)); end
            [~,slot]=min(vecnorm(elite-X(i,:),2,2)); direction=elite(slot,:)-X(i,:);
            velocity(i,:)=rand(1,D).*direction+TT*randn(1,D); TT=TT*.99;
        end
        other=randperm(N); other=other(find(other~=i,1));
        if isfinite(loss(i))&&isfinite(loss(other))&&norm(X(i,:)-X(other,:))<sigma*(maxGeneration-batch)/maxGeneration
            Px(i)=Px(i)+Pa*Px(i)*(1-Px(i)-(loss(other)/loss(i))*Px(other));
            velocity(i,:)=velocity(i,:)+((loss(other)-loss(i))/loss(other))*(X(other,:)-X(i,:));
        end
        assert(all(isfinite(velocity(i,:)))&&isfinite(Px(i)));
        if Px(i)<=.001, Px(i)=1/N; X(i,:)=rand(1,D); velocity(i,:)=0; loss(i)=Inf; end
    end
    batch=batch+1; queue=[];
    if e.bestF>before+1e-12, stagnant=0; else, stagnant=stagnant+1; end
    if enableExchange&&stagnant>=c.stagnation&&batch-lastTrigger>=c.stagnation&&e.remaining()
        [~,order]=sort(loss,'descend'); queue=order(1:max(1,ceil(c.recipientFraction*N))); lastTrigger=batch;
    end
end
result=e.finish('EPPE',seed,struct('fixedPriority',true,'donorPreference',false,'useAdd',useAdd,'useRemove',useRemove,'exchangeCount',exchangeCount));
end
function x=exchange(x,p,c,stream,e,useAdd,useRemove)
S=decode_layout(x,p); T=S; removed=[]; added=[];
request=min([1+(rand(stream)<c.twoSwapProbability) p.k numel(p.available)-p.k]);
w=.5/numel(p.ids)+.5*p.demandWeights;
for k=1:request
    outgoing=setdiff(S,removed); incoming=setdiff(p.available,union(S,added));
    if isempty(outgoing)||isempty(incoming), break; end
    preferred=intersect(p.P0,incoming); if isempty(preferred), preferred=incoming; end
    U=covered(T,p);
    if ~useAdd||rand(stream)<c.exploreProbability, a=incoming(randi(stream,numel(incoming)));
    else
        gain=sum((p.coverage(preferred,:)&~U).*w,2); [~,j]=max(gain); a=preferred(j); e.proxyCalls=e.proxyCalls+numel(preferred);
    end
    if ~useRemove, b=outgoing(randi(stream,numel(outgoing)));
    else
        uniqueWeight=zeros(size(outgoing));
        for j=1:numel(outgoing), without=covered(setdiff(T,outgoing(j)),p); uniqueWeight(j)=sum(w(U&~without)); end
        e.proxyCalls=e.proxyCalls+numel(outgoing); [~,j]=min(uniqueWeight); b=outgoing(j);
    end
    removed(end+1)=b; added(end+1)=a; T=union(setdiff(T,b),a); %#ok<AGROW>
end
for k=1:numel(added)
    i=find(p.available==removed(k)); j=find(p.available==added(k));
    if x(i)==x(j), [~,ord]=sort(x); x(ord)=((1:numel(x))-.5)/numel(x); end
    x([i j])=x([j i]);
end
assert(isequal(decode_layout(x,p),sort(T)));
end
function U=covered(S,p), U=any(p.coverage(sort([p.existing S]),:),1); end
