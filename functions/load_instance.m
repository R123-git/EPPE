function p=load_instance(c,network,condition)
a=load(fullfile(c.root,'data','instances',[network '_' condition '.mat']),'p');p=a.p;p.objective='QFloor';
assert(p.k==p.B-numel(p.existing));
assert(all(ismember(p.initialSet,p.available))&&numel(p.initialSet)==p.k);
end
