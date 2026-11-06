function [p,audit]=build_hydraulic_matrix(base,scale,c)
timer=tic; root=c.root; net=base.network; model=fullfile(root,'data','models',[net '.inp']); if strcmp(net,'Ctown'), model=fullfile(root,'data','models','CTOWN.INP'); end
work=tempname(fullfile(root,'results','hydraulics')); mkdir(work); copyfile(model,fullfile(work,'network.inp'));
old=pwd; cd(work); restore=onCleanup(@()cd(old)); %#ok<NASGU>
evalc('d=epanet(''network.inp'');'); cleanup=onCleanup(@()quietUnload(d)); %#ok<NASGU>
d.setTimeSimulationDuration(0); d.setDemandModel('DDA',0,0,0); d.setOptionsAccuracyValue(c.accuracy); d.setOptionsMaxTrials(c.maxTrials);
junc=d.getNodeJunctionIndex; [ids,ord]=sort(string(d.getNodeNameID(junc))'); junc=junc(ord); assert(isequal(ids(:),base.ids(:)));
if strcmpi(string(d.getFlowUnits),'GPM'), pf=.7032496149; ff=.0630901964; else, assert(strcmpi(string(d.getFlowUnits),'LPS')); pf=1; ff=1; end
original=d.getNodeBaseDemands; counts=d.getNodeDemandCategoriesNumber; idx=[]; cats=[]; vals=[];
for j=1:numel(junc), for cat=1:counts(junc(j))
    value=double(original{cat}(junc(j))); if value>0, idx(end+1)=junc(j); cats(end+1)=cat; vals(end+1)=value; end %#ok<AGROW>
end, end
for j=1:numel(vals), err=d.apiENsetbasedemand(idx(j),cats(j),vals(j)*scale,d.LibEPANET,d.ph); assert(err==0); end
d.addPattern('PPE2_CONST',1); multiplier=d.getOptionsPatternDemandMultiplier; assert(multiplier>0); N=numel(junc); impulseCat=zeros(N,1);
for j=1:N, impulseCat(j)=d.addNodeJunctionDemand(junc(j),0,'PPE2_CONST','PPE2_IMPULSE'); end
[pressure,demand,nominalCode]=solve(d,junc,pf,ff); assert(nominalCode==0&&all(isfinite(pressure))&&all(pressure(base.eligible)>0));
delta=nan(N,N); valid=false(N,1); errors=zeros(N,1); minimumPressure=nan(N,1); flowError=nan(N,1); impulse=base.impulseLps;
for ell=1:N
    err=d.apiENsetbasedemand(junc(ell),impulseCat(ell),impulse/ff/multiplier,d.LibEPANET,d.ph); assert(err==0);
    [changed,changedDemand,code]=solve(d,junc,pf,ff); dp=pressure-changed;
    errors(ell)=code; delta(ell,:)=dp; minimumPressure(ell)=min(changed(base.eligible)); flowError(ell)=abs(changedDemand(ell)-demand(ell)-impulse);
    valid(ell)=code==0&&all(isfinite(changed))&&all(changed(base.eligible)>0)&&abs(dp(ell))>c.pressureTolerance;
    err=d.apiENsetbasedemand(junc(ell),impulseCat(ell),0,d.LibEPANET,d.ph); assert(err==0);
    if mod(ell,50)==0||ell==N, fprintf('EXTERNAL_MATRIX %s x%.1f %d/%d invalid=%d\n',net,scale,ell,N,sum(~valid(1:ell))); end
end
[reset,~,resetCode]=solve(d,junc,pf,ff); resetScaled=max(abs(reset-pressure));
for j=1:numel(vals), err=d.apiENsetbasedemand(idx(j),cats(j),vals(j),d.LibEPANET,d.ph); assert(err==0); end
[resetOriginal,~,originalCode]=solve(d,junc,pf,ff); resetNominal=max(abs(resetOriginal-base.nominalPressure(:)'));
L=nan(N,N); diagonal=diag(delta); L(valid,:)=delta(valid,:)./diagonal(valid); Z=nan(N,N);
for col=1:N
    x=L(valid,col); sd=std(x); if sd<=1e-12, zz=zeros(size(x)); else, zz=abs((x-mean(x))/sd); end
    width=max(zz)-min(zz); if width<=1e-12, Z(valid,col)=0; else, Z(valid,col)=(zz-min(zz))/width; end
end
R=zeros(N,N); ii=find(valid);
for i=reshape(ii,1,[]), R(i,ii)=1-sqrt(mean((Z(i,:)-Z(ii,:)).^2,2))'; end
R=max(0,min(1,R)); H=abs(pressure(:)-pressure(:)'); threshold=base.settings.threshold;
C=base.distance<=threshold(1)&H<=threshold(2)&R>=threshold(3); C(~valid,:)=false; C(:,~valid)=false;
p=base; p.coverage=C; p.nominalPressure=pressure(:); p.nominalDemand=demand(:); w=max(demand,0); p.demandWeights=w/sum(w); p.deltaPressure=delta; p.R=R; p.H=H; p.validResponse=valid; p.hydraulicErrors=errors; p.flowError=flowError;
audit=struct('scale',scale,'modelHash',study_hash(model,true),'builderHash',study_hash(mfilename('fullpath')+".m",true),'nominalCode',nominalCode,'resetCode',resetCode,'originalCode',originalCode,'resetScaled',resetScaled,'resetNominal',resetNominal,'positiveCategories',numel(vals),'invalid',sum(~valid),'responseCount',N,'solveCount',N+3,'maxFlowError',max(flowError),'minimumPressure',minimumPressure,'fixedImpulseLps',impulse,'seconds',toc(timer),'builtAt',char(datetime('now','TimeZone','Asia/Shanghai')));
p.signature=study_hash(struct('base',base.signature,'audit',audit,'coverage',C));
assert(resetCode==0&&originalCode==0&&max(flowError)<max(1e-4,1e-5*impulse));
end
function [p,q,code]=solve(d,junc,pf,ff)
d.openHydraulicAnalysis; d.initializeHydraulicAnalysis(10); evalc('d.runHydraulicAnalysis;'); code=d.Errcode;
p=double(d.getNodePressure); p=reshape(p(junc)*pf,1,[]); q=double(d.getNodeActualDemand); q=reshape(q(junc)*ff,1,[]); d.closeHydraulicAnalysis;
end
function quietUnload(d), evalc('d.unload;'); end
