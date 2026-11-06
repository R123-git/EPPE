classdef RobustEvaluator < CoverageEvaluator
    methods
        function obj=RobustEvaluator(p,cap), obj@CoverageEvaluator(p,cap); end
        function f=rescore(obj,S)
            [Q,F]=robust_score(S,obj.problem); p=obj.problem; f=Q;
            if strcmp(p.objective,'F'), f=F;
            elseif strcmp(p.objective,'QFloor') && F < p.floor-p.floorTolerance
                f=F/p.floor-1;
            end
        end
        function loss=evaluate(obj,x)
            assert(obj.count<obj.cap); S=decode_layout(x,obj.problem); key=sprintf('%d_',S);
            if isKey(obj.cache,key), f=obj.cache(key); obj.cacheHits=obj.cacheHits+1;
            else, f=obj.rescore(S); obj.cache(key)=f; obj.actualScores=obj.actualScores+1; end
            if f>obj.bestF+1e-12, obj.bestF=f; obj.bestKeys=x; obj.bestSet=S; end
            obj.count=obj.count+1; obj.trace(obj.count)=obj.bestF; obj.time(obj.count)=toc(obj.startClock); loss=2-f;
        end
        function [children,losses]=exchangeBatch(obj,x,S,removed,incoming)
            incoming=incoming(1:min(numel(incoming),obj.cap-obj.count)); children=repmat(x,numel(incoming),1); losses=zeros(numel(incoming),1);
            for j=1:numel(incoming)
                i=find(obj.problem.available==removed); k=find(obj.problem.available==incoming(j));
                if children(j,i)==children(j,k), [~,ord]=sort(children(j,:)); children(j,ord)=((1:numel(x))-.5)/numel(x); end
                children(j,[i k])=children(j,[k i]); assert(isequal(decode_layout(children(j,:),obj.problem),sort([setdiff(S,removed) incoming(j)])));
                losses(j)=obj.evaluate(children(j,:));
            end
        end
        function r=finish(obj,method,seed,mechanisms)
            p=obj.problem; S=obj.bestSet; [Q,F,loss,worst]=robust_score(S,p); target=obj.rescore(S);
            assert(abs(target-obj.bestF)<1e-12&&obj.actualScores+obj.cacheHits==obj.count&&numel(unique(S))==p.k&&all(diff(obj.trace(1:obj.count))>=-1e-12));
            total=sort([p.existing S]); r=struct('method',method,'seed',seed,'target',target,'Q',Q,'F',F,'worstLoss',loss,'worstSensor',p.ids(total(worst)),'S',S,'queries',obj.count,'actualScores',obj.actualScores,'cacheHits',obj.cacheHits,'proxyCalls',obj.proxyCalls,'seconds',toc(obj.startClock),'trace',obj.trace(1:obj.count),'mechanisms',mechanisms);
            r.feasible=F>=p.floor-p.floorTolerance; r.floor=p.floor; r.Fstar=p.Fstar; r.initialSet=p.initialSet;
        end
    end
end
