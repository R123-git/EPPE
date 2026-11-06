function net=read_network(file)
lines=splitlines(string(fileread(file))); section=""; ids=strings(0); xy=zeros(0,2); links=cell(0,3); vertices=cell(0,3); reservoirs=strings(0); tanks=strings(0);
for line=lines'
    line=strtrim(extractBefore(line+";",";")); if strlength(line)==0, continue; end
    if startsWith(line,"["), section=upper(line); continue; end
    tok=split(line); tok(tok=="")=[];
    switch section
        case "[COORDINATES]", ids(end+1)=tok(1); xy(end+1,:)=[str2double(tok(2)) str2double(tok(3))];
        case {"[PIPES]","[PUMPS]","[VALVES]"}, links(end+1,:)={tok(1),tok(2),tok(3)};
        case "[VERTICES]", vertices(end+1,:)={tok(1),str2double(tok(2)),str2double(tok(3))};
        case "[RESERVOIRS]", reservoirs(end+1)=tok(1);
        case "[TANKS]", tanks(end+1)=tok(1);
    end
end
assert(all(isfinite(xy),'all')&&numel(unique(ids))==numel(ids));
xx=[]; yy=[];
for j=1:size(links,1)
    [ok,ix]=ismember([string(links{j,2}) string(links{j,3})],ids); assert(all(ok));
    v=find(string(vertices(:,1))==string(links{j,1})); mid=zeros(0,2); if ~isempty(v), mid=cell2mat(vertices(v,2:3)); end
    pts=[xy(ix(1),:);mid;xy(ix(2),:)]; xx=[xx;pts(:,1);NaN]; yy=[yy;pts(:,2);NaN]; %#ok<AGROW>
end
net=struct('ids',ids,'xy',xy,'xline',xx,'yline',yy,'reservoirs',reservoirs,'tanks',tanks,'links',size(links,1));
end
