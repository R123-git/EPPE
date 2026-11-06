function records=load_records(c,source)
if strcmp(source,'archive')
 a=load(fullfile(c.root,'data','archive','records.mat'),'records');records=a.records;
elseif strcmp(source,'new')
 files=dir(fullfile(c.root,'results','runs','*.mat'));records=cell(1,numel(files));
 for j=1:numel(files),a=load(fullfile(files(j).folder,files(j).name),'r');records{j}=a.r;end
else,error('Source must be archive or new.');end
end
