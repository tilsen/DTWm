function [] = fig_process_localrate()

dbstop if error; close all
h = tw_helpers;

ax = stf([3 1],[0.11 0.065 0.01 0.025],[0 0.025],'aspect',0.85);

S = test_signals('local_rate_example');

dt = S.dt(1);
t = S.t{1};
y = S.X{1};
theta = S.state{1};

colors = [0 0 0; 1 0 0];
ms = 4;

for i=1:2

    S.est_state{i} = -log(S.X{i})/2;
    S.est_dstate{i} = diff(S.est_state{i});
    S.est_rate{i} = S.est_dstate{i}/dt;

    axes(ax(1));
    plot(S.t{i},S.X{i},'o-','linew',2,'markersize',ms,'color',colors(i,:),'markerfacecolor','w'); hold on;

    axes(ax(2));
    plot(S.t{i},S.est_state{i},'o-','linew',2,'markersize',ms,'color',colors(i,:),'markerfacecolor','w'); hold on;

    axes(ax(3));
    plot(S.t{i}(1:end-1),S.est_rate{i},'o-','linew',2,'markersize',ms,'color',colors(i,:),'markerfacecolor','w'); hold on;    
end


%%
ix = find(t>=0.15,1,'first');
dstate = S.est_dstate{1}(ix);
rate = S.est_rate{1}(ix);

axis(ax,'tight');
axrescalex(0.01,ax);
axrescaley(0.05,ax);
set(ax,'YGrid','on','XGrid','on','Box','off','ticklen',0.003*[1 1],'fontsize',h.fs(end));
set(ax(1:end-1),'XTickLabel',[]);

%%

xof = @(th)sum(get(th,'Extent').*[1 0 1 0]);
xot = 0.05;
yo = 0.025;

axes(ax(1));

bh(1) = drawbrace([t(ix) y(ix)+yo],[t(ix+1) y(ix)+yo],0.0035,'color','k','linew',1);

str = ['$\Delta t=' num2str(dt,'%1.3f') '\ $'];
th(1) = text(t(ix),max(bh(1).YData),str,'hori','left','verti','bot','fontsize',h.fs(3),'interp','latex');

str = 'process output:';
th(end+1) = text(xot,max(ylim),str,'hori','left','verti','top','fontsize',h.fs(2));

str = '$\ \ {\it{y}}(t) = e^{-2\theta(t)}$';
th(end+1) = text(xof(th(end)),max(ylim),str,'hori','left','verti','top','fontsize',h.fs(2),'interp','latex');

ixo = ix+4;
th(end+1) = text(t(ixo),y(ixo)-0.05,'variable-rate process','hori','right','verti','top','fontsize',h.fs(2),'color','r');
th(end+1) = text(t(ixo),y(ixo),'constant-rate process','hori','left','verti','bot','fontsize',h.fs(2));


%%
axes(ax(2));

bh(2) = drawbrace([t(ix+1) theta(ix)-yo],[t(ix) theta(ix)-yo],0.0035,'color','k','linew',1);

str = 'process state:';
th(end+1) = text(xot,max(ylim),str,'hori','left','verti','top','fontsize',h.fs(2));

str = '$\ \ \theta(t) = \theta(t-\Delta t)+r(t)\Delta t$';
th(end+1) = text(xof(th(end)),max(ylim),str,'hori','left','verti','top','fontsize',h.fs(2),'interp','latex');

str = ['$\Delta \tilde{\theta}(t) = \frac{1}{2}(-\log{y(t+\Delta t)}+\log{y(t)}) = ' num2str(dstate,'%1.3f') '$'];
th(end+1) = text(t(ix)-0.025, min(bh(2).YData)-0.025,str,'FontSize',h.fs(2)-2,'interp','latex','verti','top');

%%

axes(ax(3));

%bh(3) = drawbrace([t(ix+1) S.est_rate{1}(ix)-yo],[t(ix) S.est_rate{1}(ix)-yo],0.0035,'color','k','linew',1);

%str = ['$\tilde{r}(t)=' num2str(rate,'%1.3f') '\ $'];
%th(end+1) = text(t(ix),max(bh(3).YData)-yo,str,'hori','left','verti','top','fontsize',h.fs(3),'interp','latex');

str = 'inferred local rate:';
th(end+1) = text(xot,max(ylim),str,'hori','left','verti','top','fontsize',h.fs(2));

str = ['$\ \ \tilde{r}(t)=\frac{\Delta \tilde{\theta}(t)}{\Delta t}$'];
th(end+1) = text(xof(th(end))+0.01,max(ylim),str,'FontSize',h.fs(2),'interp','latex','verti','top');

ax(3).XLim = ax(2).XLim;

%%

xlabel(ax(end),'time (t)','fontsize',h.fs(2));
ylabel(ax(1),'y(t) ','fontsize',h.fs(2),'hori','right','rotation',0);
ylabel(ax(2),'\theta(t) ','fontsize',h.fs(2),'hori','right','rotation',0);
ylabel(ax(3),'r(t) ','fontsize',h.fs(2),'hori','right','rotation',0);

stfig_panlab(ax,{'A','B','C'},'fontsize',h.fs(1),'verti','mid','hori','right','xoffset',-0.075);

%%
h.printfig(mfilename);

end




