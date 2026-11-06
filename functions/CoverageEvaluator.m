classdef CoverageEvaluator < handle
    properties
        problem; cap; count=0; actualScores=0; cacheHits=0;
        bestF=-Inf; bestKeys=[]; bestSet=[]; trace; time;
        cache; startClock; proxyCalls=0;
    end
    methods
        function obj=CoverageEvaluator(p,cap)
            obj.problem=p; obj.cap=cap;
            obj.trace=zeros(cap,1); obj.time=zeros(cap,1);
            obj.cache=containers.Map('KeyType','char','ValueType','double');
            obj.startClock=tic;
        end
        function loss=evaluate(obj,x)
            assert(obj.count<obj.cap,'FE预算已耗尽');
            S=decode_layout(x,obj.problem); key=sprintf('%d_',S);
            if isKey(obj.cache,key)
                F=obj.cache(key); obj.cacheHits=obj.cacheHits+1;
            else
                F=coverage_score(S,obj.problem); obj.cache(key)=F;
                obj.actualScores=obj.actualScores+1;
            end
            if F>obj.bestF+1e-12
                obj.bestF=F; obj.bestKeys=x; obj.bestSet=S;
            end
            obj.count=obj.count+1; obj.trace(obj.count)=obj.bestF;
            obj.time(obj.count)=toc(obj.startClock); loss=2-F;
        end
        function yes=remaining(obj), yes=obj.count<obj.cap; end
        function r=finish(obj,method,seed,mechanisms)
            [F,CN,CW]=coverage_score(obj.bestSet,obj.problem);
            assert(abs(F-obj.bestF)<1e-12 && numel(obj.bestSet)==obj.problem.k);
            r=struct('method',method,'seed',seed,'F',F,'CN',CN,'CW',CW, ...
                'newNodes',obj.problem.ids(obj.bestSet),'existingNodes',obj.problem.ids(obj.problem.existing), ...
                'newIndices',obj.bestSet,'totalBudget',obj.problem.B,'queries',obj.count, ...
                'actualScores',obj.actualScores,'cacheHits',obj.cacheHits, ...
                'proxyCalls',obj.proxyCalls,'seconds',toc(obj.startClock), ...
                'trace',obj.trace(1:obj.count),'time',obj.time(1:obj.count), ...
                'mechanisms',mechanisms,'instanceHash',obj.problem.signature);
        end
    end
end
