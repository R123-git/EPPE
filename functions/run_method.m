function r=run_method(method,p,c,seed)
switch method
 case {'EPPE','PPE','PPE_M1','PPE_M2','RandomExchange'}
  r=EPPE(p,c,seed,method);
 case 'EROA',r=alg_eroa(p,c,seed);
 case {'SSA','SOA','DBO'},r=alg_reference(p,c,seed,method);
 otherwise,error('Unknown method: %s',method);
end
r.method=method;
end
