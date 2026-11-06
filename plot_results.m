function plot_results(source)
% Four manuscript figures. Existing data are never smoothed or truncated.
if nargin<1,source='archive';end
c=config();records=load_records(c,source);out=fullfile(c.root,'results',source,'figures');
if ~isfolder(out),mkdir(out);end
cases={'Ctown','H1';'Modena','M1'};
colors=[.85 .1 .12;0 .35 .75;.1 .6 .3;.7 .2 .7;0 .65 .75;.18 .18 .18];
styles={'-','--','-.',':','--','-.'};marks={'o','s','^','d','v','>'};
for kind=1:4
 f=figure('Color','w','Visible','off','Units','centimeters','Position',[2 2 18 25]);
 tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');
 for k=1:2
  net=cases{k,1};cond=cases{k,2};p=load_instance(c,net,cond);ax=nexttile;hold(ax,'on');
  label=net;if strcmp(net,'Ctown'),label='C-Town';end
  if kind==1||kind==3
   names=c.main;display=names;cc=colors;st=styles;mk=marks;
   if kind==3,names={'PPE','PPE_M1','PPE_M2','EPPE'};display={'PPE','PPE + M1','PPE + M2','EPPE'};cc=colors([2 3 4 1],:);st={'--','-.',':','-'};mk={'o','s','^','d'};end
   for m=1:numel(names)
    rr=select(records,net,cond,names{m});tr=cellfun(@(r)r.trace(:),rr,'UniformOutput',false);tr=100*cat(2,tr{:});
    assert(all(tr>=0,'all')&&all(diff(tr)>=-1e-10,'all'));
    plot(ax,1:size(tr,1),mean(tr,2),'Color',cc(m,:),'LineStyle',st{m},'Marker',mk{m},'MarkerIndices',unique([1 500:500:size(tr,1)]),'MarkerSize',3,'LineWidth',1.2,'DisplayName',display{m});
   end
   xlabel(ax,'Complete objective queries (FE)');ylabel(ax,'Mean best feasible Q (%)');xlim(ax,[1 4000]);
   if kind==1,ylim(ax,[80 100]);end
   legend(ax,'Location','northwest','FontSize',8,'Interpreter','none');clear_corner(ax);
  elseif kind==2
   for m=1:6
    rr=select(records,net,cond,c.main{m});q=100*cellfun(@(r)r.Q,rr);q=q(:);
    boxchart(ax,repmat(m,numel(q),1),q,'BoxFaceColor',colors(m,:),'MarkerStyle','none','HandleVisibility','off');
    scatter(ax,m+linspace(-.1,.1,numel(q))',q,10,colors(m,:),'filled','MarkerFaceAlpha',.55,'HandleVisibility','off');
   end
   xticks(ax,1:6);xticklabels(ax,c.main);xlim(ax,[.5 6.5]);ylabel(ax,'Final worst-failure Q (%)');
   yl=ylim(ax);ylim(ax,[yl(1) yl(2)+.2*diff(yl)]);
   plot(ax,NaN,NaN,'ks','DisplayName','Box: IQR and median');plot(ax,NaN,NaN,'k.','DisplayName','Points: individual runs');
   legend(ax,'Location','northeast','FontSize',8);
  else
   rr=select(records,net,cond,'EPPE');q=cellfun(@(r)r.Q,rr);runs=cellfun(@(r)r.run,rr);q=q(:);runs=runs(:);
   [~,ord]=sortrows([q runs],[1 2]);r=rr{ord(max(1,floor(numel(rr)/2)))};
   [v,~,T]=layout_metrics(r.S,p);assert(max(abs([v.Q v.F]-[r.Q r.F]))<1e-10);
   model='Modena.inp';if strcmp(net,'Ctown'),model='CTOWN.INP';end
   g=read_network(fullfile(c.root,'data','models',model));[ok,pos]=ismember(p.ids,g.ids);assert(all(ok));
   origin=min(g.xy,[],1);g.xy=g.xy-origin;g.xline=g.xline-origin(1);g.yline=g.yline-origin(2);
   plot(ax,g.xline,g.yline,'-','Color',[.7 .76 .84],'LineWidth',.4,'HandleVisibility','off');
   scatter(ax,g.xy(pos,1),g.xy(pos,2),3,[.65 .72 .81],'filled','HandleVisibility','off');
   counts=sum(p.coverage(T,:),1);[~,critical]=ismember(v.worstSensor,p.ids);
   lost=find(p.coverage(critical,:)&counts==1);backup=find(counts>=2);
   scatter(ax,g.xy(pos(T),1),g.xy(pos(T),2),21,'r','filled','DisplayName','Installed sensors');
   scatter(ax,g.xy(pos(p.existing),1),g.xy(pos(p.existing),2),39,'ks','LineWidth',1,'DisplayName','Existing sensors');
   xy=g.xy(pos(critical),:);scatter(ax,xy(1),xy(2),100,'p','MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','Critical sensor');
   scatter(ax,g.xy(pos(lost),1),g.xy(pos(lost),2),22,[.9 .35 0],'x','LineWidth',1,'DisplayName','Lost targets');
   scatter(ax,g.xy(pos(backup),1),g.xy(pos(backup),2),18,[0 .5 .6],'o','DisplayName','Backup targets (>=2)');
   for type=1:2
    list=g.reservoirs;mark='s';clr=[.2 .65 .3];if type==2,list=g.tanks;mark='^';clr=[0 .65 .8];end
    ix=find(ismember(g.ids,list));scatter(ax,g.xy(ix,1),g.xy(ix,2),32,clr,mark,'filled','HandleVisibility','off');
   end
   tag='new';if ismember(critical,p.existing),tag='existing';end
   lo=min(g.xy,[],1);span=max(g.xy,[],1)-lo;
   labelxy=lo+[.02 .16].*span;location='northwest';
   if strcmp(net,'Modena'),labelxy=lo+[.02 .91].*span;location='northeast';end
   plot(ax,[labelxy(1)+.25*span(1) xy(1)],[labelxy(2) xy(2)],'k-','LineWidth',.6,'HandleVisibility','off');
   text(ax,labelxy(1),labelxy(2),sprintf('%s (%s)\n%.3f pp; %d targets',char(v.worstSensor),tag,100*(r.F-r.Q),numel(lost)),'FontSize',8,'VerticalAlignment','bottom','BackgroundColor','w','Margin',2);
   axis(ax,'equal');axis(ax,'padded');xlabel(ax,'Relative model X coordinate');ylabel(ax,'Relative model Y coordinate');
   legend(ax,'Location',location,'FontSize',7);
  end
  title(ax,sprintf('(%c) %s (%s)','a'+k-1,label,cond));grid(ax,'on');box(ax,'on');
  set(ax,'FontName','Times New Roman','FontSize',9);
 end
 names={'Algorithm convergence','Final coverage distributions','Component convergence','Sensor layouts and critical outages'};
 set(f,'Name',names{kind},'NumberTitle','off','InvertHardcopy','off');
 stem=fullfile(out,sprintf('Figure_%d',kind));
 exportgraphics(f,[stem '.tif'],'Resolution',600,'Colorspace','rgb');savefig(f,[stem '.fig']);close(f);
end
fprintf('Saved Figures 1-4 under results/%s/figures.\n',source);
end

function rr=select(records,net,cond,method)
use=cellfun(@(r)strcmp(r.network,net)&&strcmp(r.condition,cond)&&strcmp(r.method,method),records);rr=records(use);
assert(~isempty(rr),'Missing results: %s %s %s',net,cond,method);
[~,order]=sort(cellfun(@(r)r.run,rr));rr=rr(order);
end

function clear_corner(ax)
% Select an internal legend corner using plotted segments and target markers.
lg=legend(ax);lg.Units='normalized';ax.Units='normalized';locations={'northwest','southeast','northeast','southwest'};
for attempt=1:5
 best=Inf;chosen='northwest';
 for j=1:4
  lg.Location=locations{j};drawnow;ap=ax.InnerPosition;lp=lg.Position;
  rect=[(lp(1:2)-ap(1:2))./ap(3:4) lp(3:4)./ap(3:4)];xy=[];
  for h=allchild(ax)'
   if ~isprop(h,'XData')||~isprop(h,'YData'),continue;end
   x=h.XData(:);y=h.YData(:);if numel(x)~=numel(y),continue;end
   x=(x-ax.XLim(1))/diff(ax.XLim);y=(y-ax.YLim(1))/diff(ax.YLim);
   if isprop(h,'LineStyle')&&~strcmp(h.LineStyle,'none')
    xx=[];yy=[];for k=1:numel(x)-1,xx=[xx;linspace(x(k),x(k+1),15)'];yy=[yy;linspace(y(k),y(k+1),15)'];end %#ok<AGROW>
    x=[x;xx];y=[y;yy];
   end
   xy=[xy;x y]; %#ok<AGROW>
  end
  hits=sum(xy(:,1)>=rect(1)-.012&xy(:,1)<=rect(1)+rect(3)+.012&xy(:,2)>=rect(2)-.012&xy(:,2)<=rect(2)+rect(4)+.012);
  if hits<best,best=hits;chosen=locations{j};end
 end
 lg.Location=chosen;if best==0,return;end
 yl=ax.YLim;ax.YLim=[yl(1) yl(2)+.15*diff(yl)];
end
end
