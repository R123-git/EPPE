function r=alg_eroa(p,cfg,seed)
rng(seed,'twister'); e=RobustEvaluator(p,cfg.FE.(p.network));
N=cfg.population; D=numel(p.available); primesList=primes(100*D);
prime=primesList(find(primesList>=2*D+3,1));
X=-1+2*mod((1:N)'*(2*cos(2*pi*(1:D)/prime)),1);
X(1,:)=2*p.initialKeys-1; leaders=zeros(0,D); leaderLoss=[]; v2=2*pi*rand(N,D); removed=false(1,10);
while e.remaining()
    loss=Inf(N,1);
    for i=1:N
        if ~e.remaining(), break; end
        loss(i)=e.evaluate((X(i,:)+1)/2);
    end
    allX=[leaders;X]; allLoss=[leaderLoss;loss]; [allLoss,ord]=sort(allLoss);
    take=min(3,sum(isfinite(allLoss))); leaders=allX(ord(1:take),:); leaderLoss=allLoss(1:take);
    while size(leaders,1)<3, leaders(end+1,:)=leaders(1,:); leaderLoss(end+1)=leaderLoss(1); end
    pool=[leaders;mean(leaders,1)];
    if ~e.remaining(), break; end
    u=e.count/e.cap; cycle=min(9,floor(10*u)); phase=10*u-cycle;
    [~,order]=sort(loss);
    if u>=(cycle+.5)/10 && ~removed(cycle+1)
        worst=order(end); X(worst,:)=[]; v2(worst,:)=[]; loss(worst)=[];
        N=N-1; removed(cycle+1)=true; [~,order]=sort(loss);
    end
    if phase<.5
        for rank=1:N
            i=order(rank); target=pool(randi(4),:);
            if rank<=N/3
                a=randi(N); radius=norm(X(a,:)-target);
                X(i,:)=target+radius*cos(pi*rand(1,D)).*sin(pi/2*rand(1,D));
            else
                C=2*rand-1; theta=2*pi*rand(1,D); w=1/40;
                v1=sqrt(2.5^2*w^2*sin(w+theta).^2+.1^2*w^2*cos(w+theta).^2)/2;
                v2(i,:)=v2(i,:).*w.*cos(w+theta-.78545);
                X(i,:)=X(i,:)+C*(v1+v2(i,:))+(target-X(i,:))*(1-C).*rand(1,D);
            end
        end
    else
        perm=randperm(3); groupSize=round(N/3);
        for i=1:N
            group=min(3,ceil(i/groupSize)); target=pool(perm(group),:);
            if group==1, step=(2*rand(1,D)-1)*exp(u-1).*sign(rand(1,D)-.5);
            elseif group==2
                beta=1.5; sd=(gamma(1+beta)*sin(pi*beta/2)/(beta*gamma((1+beta)/2)*2^((beta-1)/2)))^(1/beta);
                step=sd*randn(1,D)./max(abs(randn(1,D)).^(1/beta),realmin);
            else, step=.1*tan(pi*(rand(1,D)-.5)); end
            X(i,:)=target+step;
        end
    end
    X=max(-1,min(1,X));
end
r=e.finish('EROA',seed,struct('budgetAdapted',true));
end
