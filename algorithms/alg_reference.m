function r=alg_reference(p,c,seed,name)
rng(seed,'twister'); e=RobustEvaluator(p,c.FE.(p.network)); N=c.population; D=numel(p.available);
if strcmp(name,'SSA'), M=ceil((e.cap-N)/(N+20)); elseif strcmp(name,'DBO'), M=ceil((e.cap-N)/N); else, M=ceil(e.cap/N); end
try
    switch name
        case 'SSA', ssa_engine(N,M,0,1,D,@evaluate,p.initialKeys);
        case 'SOA', soa_paper_engine(N,M,0,1,D,@evaluate,p.initialKeys);
        case 'DBO', dbo_engine(N,M,0,1,D,@evaluate,p.initialKeys);
        otherwise, error('Unknown reference');
    end
catch ex
    if ~strcmp(ex.identifier,'EPPE:Budget'), rethrow(ex); end
end
assert(e.count==e.cap); r=e.finish(name,seed,struct('downloadedKernel',~strcmp(name,'SOA'),'paperEquations',strcmp(name,'SOA'),'maximumIterations',M,'sharedWarmStart',true));
    function loss=evaluate(x)
        if ~e.remaining(), error('EPPE:Budget','Shared FE cap reached'); end
        loss=e.evaluate(x);
    end
end
