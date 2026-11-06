function rebuilt=prepare_data(toolkitFolder)
% Optional hydraulic rebuild. Frozen placement instances are never overwritten.
c=config();if nargin>0,addpath(genpath(toolkitFolder));end
assert(exist('epanet','file')==2,'Provide the EPANET MATLAB toolkit folder.');
if ~isfolder(c.dataDir),mkdir(c.dataDir);end
c.toolkit=fileparts(which('epanet'));rebuilt=cell(1,2);
for j=1:2
 rebuilt{j}=make_problem(c.networks{j},c);
 for ci=1:2
  p=make_instance(rebuilt{j},sprintf('R%d',ci),c);p.condition=c.conditions{j}{ci};
  [p,certificate]=normal_solution(p,c.B.(p.network));
  save(fullfile(c.dataDir,[p.network '_' p.condition '.mat']),'p','certificate');
 end
end
fprintf('Rebuilt hydraulic matrices saved under results/hydraulics.\n');
end
