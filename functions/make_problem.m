function p=make_problem(network,cfg)
if ~isfolder(cfg.dataDir), mkdir(cfg.dataDir); end
assert(exist('epanet','file')==2,'Add the EPANET MATLAB toolkit to the path first.');
file=fullfile(cfg.root,'data','models',[network '.inp']); if strcmp(network,'Ctown'),file=fullfile(cfg.root,'data','models','CTOWN.INP');end
assert(isfile(file),'找不到模型：%s',file);
settings=struct('model',study_hash(file,true),'builder',study_hash(mfilename('fullpath')+".m",true), ...
    'threshold',cfg.threshold.(network),'accuracy',cfg.accuracy,'trials',cfg.maxTrials, ...
    'impulse',cfg.impulseFraction,'tolerance',cfg.pressureTolerance,'invalid',cfg.maxInvalidFraction, ...
    'epanetClass',study_hash(which('epanet'),true));
dll=dir(fullfile(fileparts(which('epanet')),'64bit','*.dll'));
settings.dllHashes=cell(numel(dll),1);
for k=1:numel(dll), settings.dllHashes{k}=study_hash(fullfile(dll(k).folder,dll(k).name),true); end
signature=study_hash(settings); destination=fullfile(cfg.dataDir,[network '_coverage.mat']);
if isfile(destination)
    a=load(destination,'p');
    if strcmp(a.p.signature,signature), p=a.p; return; end
end
clock=tic; nominal=simulate(file,cfg,[],0); assert(nominal.code==0,'名义水力求解未收敛');
[ids,order]=sort(string(nominal.ids(:))); junc=nominal.junc(order);
[pf,ff,lf]=units(nominal.units); pressure=nominal.pressure(junc)*pf;
demand=nominal.demand(junc)*ff; eligible=find(isfinite(pressure)&pressure>0);
N=numel(junc); assert(numel(eligible)>=cfg.B.(network));
positive=demand(demand>0); assert(~isempty(positive)); impulse=cfg.impulseFraction*median(positive);
delta=nan(N,N); valid=false(N,1); errors=zeros(N,1); flowError=nan(N,1); minimumPressure=nan(N,1);
fprintf('Building %s coverage (%d node responses)\n',network,N);
for ell=1:N
    s=simulate(file,cfg,junc(ell),impulse/ff); errors(ell)=s.code;
    changed=s.pressure(junc)*pf; dp=pressure-changed;
    flowError(ell)=abs((s.demand(junc(ell))*ff-demand(ell))-impulse);
    valid(ell)=s.code==0 && all(isfinite(changed)) && all(changed(eligible)>0) && abs(dp(ell))>cfg.pressureTolerance;
    delta(ell,:)=dp(:)'; minimumPressure(ell)=min(changed(eligible));
    if mod(ell,25)==0 || ell==N, fprintf('  %d/%d, invalid=%d\n',ell,N,sum(~valid(1:ell))); end
end
if mean(~valid)>cfg.maxInvalidFraction
    audit=struct('network',network,'ids',ids,'valid',valid,'errors',errors,'minimumPressure',minimumPressure, ...
        'response',delta,'flowError',flowError,'settings',settings,'signature',signature);
    save(fullfile(cfg.dataDir,[network '_failed_hydraulic_audit.mat']),'audit');
    error('无效扰动比例%.2f%%；停止并保存水力审计，不能静默改阈值',100*mean(~valid));
end
assert(sum(valid)>=2 && max(flowError)<max(1e-4,1e-5*impulse),'新增出流或有效样本核验失败');
L=nan(N,N);
diagonal=diag(delta); L(valid,:)=delta(valid,:)./diagonal(valid);
Z=nan(N,N);
for col=1:N
    x=L(valid,col); sd=std(x);
    if sd<=1e-12, z=zeros(size(x)); else, z=abs((x-mean(x))/sd); end
    width=max(z)-min(z);
    if width<=1e-12, Z(valid,col)=0; else, Z(valid,col)=(z-min(z))/width; end
end
R=zeros(N,N); ii=find(valid);
for i=reshape(ii,1,[]), R(i,ii)=1-sqrt(mean((Z(i,:)-Z(ii,:)).^2,2))'; end
R=max(0,min(1,R)); H=abs(pressure-pressure');
ends=nominal.ends; G=graph(ends(:,1),ends(:,2),nominal.lengths*lf,nominal.nodeCount);
D=distances(G,junc,junc); threshold=cfg.threshold.(network);
C=D<=threshold(1)&H<=threshold(2)&R>=threshold(3);
C(~valid,:)=false; C(:,~valid)=false;
weights=max(demand,0); weights=weights(:)'/sum(weights);
p=struct('network',network,'ids',ids,'coverage',C,'demandWeights',weights, ...
    'nominalPressure',pressure,'nominalDemand',demand,'eligible',eligible(:)', ...
    'distance',D,'graph',G,'junctionIndex',junc,'deltaPressure',delta, ...
    'R',R,'H',H,'validResponse',valid,'hydraulicErrors',errors,'flowError',flowError, ...
    'impulseLps',impulse,'signature',signature,'settings',settings,'buildSeconds',toc(clock));
save(destination,'p','-v7.3');
end
function s=simulate(file,cfg,target,flow)
root=fullfile(cfg.dataDir,'runtime'); if ~isfolder(root), mkdir(root); end
temporary=tempname(root); mkdir(temporary); copyfile(file,fullfile(temporary,'network.inp'));
old=pwd; cd(temporary); restore=onCleanup(@()cd(old)); %#ok<NASGU>
evalc('d=epanet(''network.inp'');'); close=onCleanup(@()quietUnload(d)); %#ok<NASGU>
d.setTimeSimulationDuration(0); d.setDemandModel('DDA',0,0,0);
d.setOptionsAccuracyValue(cfg.accuracy); d.setOptionsMaxTrials(cfg.maxTrials);
if ~isempty(target)
    d.addPattern('CCPPE_CONST',1);
    multiplier=d.getOptionsPatternDemandMultiplier; assert(multiplier>0);
    d.addNodeJunctionDemand(target,flow/multiplier,'CCPPE_CONST','CCPPE_impulse');
end
d.openHydraulicAnalysis; d.initializeHydraulicAnalysis(0); d.runHydraulicAnalysis;
s.code=d.Errcode; s.pressure=double(d.getNodePressure); s.demand=double(d.getNodeActualDemand);
s.pressure=s.pressure(:); s.demand=s.demand(:); d.closeHydraulicAnalysis;
if isempty(target)
    s.junc=d.getNodeJunctionIndex; s.ids=d.getNodeNameID(s.junc);
    s.units=d.getFlowUnits; s.ends=d.NodesConnectingLinksIndex;
    s.lengths=double(d.getLinkLength); s.lengths=s.lengths(:); s.nodeCount=d.getNodeCount;
end
end
function quietUnload(d), evalc('d.unload;'); end
function [pf,ff,lf]=units(u)
if iscell(u), u=u{1}; end
switch upper(string(u))
    case 'GPM', pf=.7032496149; ff=.0630901964; lf=.3048;
    case 'LPS', pf=1; ff=1; lf=1;
    case 'CMH', pf=1; ff=1/3.6; lf=1;
    case 'CMD', pf=1; ff=1/86.4; lf=1;
    otherwise, error('当前未定义该单位转换：%s',string(u));
end
end
