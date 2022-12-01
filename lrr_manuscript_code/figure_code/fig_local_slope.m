function [] = fig_local_slope()

dbstop if error; close all;
h = tw_helpers;

winsize = 101; %local slope window size
winc = 250; %window center for zoom

S = test_signals('growth_decay_example');

D.X = S.X([1 end]);
D.rates = S.rates([1 end]);
D.rates = cellfun(@(c){c/S.dt},D.rates);

D = repmat(D,2,1);
D(2).X = D(1).X([2 1]);
D(2).rates = D(1).rates([2 1]);
%
for i=1:length(D)
    [D(i).map,D(i).distance,D(i).map_info] = dtwm(D(i).X{1},D(i).X{2});
    [D(i).Xa(1,:),D(i).Xa(2,:),D(i).map_ref] = align_signals(D(i).map,D(i).X{1},D(i).X{2},'aligntype','reference');
    [D(i).LRR] = lrr(D(i).map,'winsize',winsize);
    D(i).relrate = D(i).rates{1}./D(i).rates{2}(D(i).map_ref(2,:));
end

map = D(1).map;
xrng = winc+winsize*[-1 1]/2;
yrng = [map(2,find(map(1,:)>=xrng(1),1,'first')) map(2,find(map(1,:)<=xrng(2),1,'last'))];

ycent(1) = round(mean(map(2,map(1,:)==winc)));
ycent(2) = round(mean(map(1,map(2,:)==ycent(1))));

winc = [winc ycent(1)]; %append window center for y

ixs_x = map(1,:)>=xrng(1) & map(1,:)<=xrng(2);
ixs_y = map(2,:)>=yrng(1) & map(2,:)<=yrng(2); 

ixs = ixs_x & ixs_y;

fprintf('i = %1.0f, matches j = %1.0f\n',winc(1),ycent(1));

%%
axpan = [1 2; 4 3];
ax = stf(axpan,[0.085 0.075 0.01 0.05],[0.10 0.125],'aspect',0.95);

colors = lines(2);
gcol = [.5 .5 .5];
lw = 2;
ls = {'-','-'};
ms = 4;

axs = stfig_subaxpos(ax(1),[1; 2],[0 0 0 0 0 0.01]);
delete(ax(1));
ax = [axs(1) ax(2:end) axs(2:end)];

%% orignal signals
axes(ax(1));
for j=1:2
    ph(j) = plot(D(j).X{1},ls{j},'color',colors(j,:),'linew',lw); hold on;
    plot(winc(j),D(j).X{1}(winc(j)),'o','linew',2,'color',colors(j,:));

end

%% process rates
axes(ax(5));
for j=1:2
    phr(j) = plot(D(j).rates{1},'-','color',colors(j,:),'linew',lw); hold on;
    plot(winc(j),D(j).rates{1}(winc(j)),'o','linew',2,'color',colors(j,:));
end

%% process relative rates
% axes(ax(6));
% for j=1:2
%     phr(j) = plot(D(j).relrate,'--','color',colors(j,:),'linew',lw); hold on;
%     plot(winc(j),D(j).relrate(winc(j)),'o','linew',2,'color',colors(j,:));
% end

%% distance matrix / warping curve
axes(ax(2));
H = plot_dtw_matrix(D(1).map_info,D(1).X{1},D(1).X{2});
set([H.s1 H.s2],'linew',2);
set(H.maph,'linew',2);

delete(H.cbh);
line(xrng([1 2 2 1 1]),yrng([1 1 2 2 1]),'color','r','linew',2);
%axis(H.axd,'equal');

%% local slope at time index i
axes(ax(3));

plot(D(1).map(1,:),D(1).map(2,:),'o','color',gcol,'markersize',ms); hold on;
plot(D(1).map(1,ixs),D(1).map(2,ixs),'ko','markersize',ms,'markerfacecolor',gcol); hold on;

loc_slope = D(1).LRR(winc(1));
loc_intercept = mean(D(1).map(2,ixs));

regfcn = @(x)loc_intercept+loc_slope*(x-winc(1));

x = unique(D(1).map(1,ixs));
y = regfcn(x);

plot(x,y,'r-','linew',3);
text(mean(x),mean(y),sprintf('  slope = %1.2f',loc_slope), ...
    'fontsize',h.fs(3),'verti','top');

xlim(xrng); ylim(yrng);
axrescale(0.25,0.25);

axis equal;
plot(winc(1)*[1 1],ylim,'k--');
plot(xlim,ycent(1)*[1 1],'k--');

xc = mean(xlim);
yc = mean(ylim);

oo = [-10 10];
arrow([xc yc] + oo,[xc yc],'tipangle',30,'length',7);

str = {['$(\phi_x(k),\phi_y(k))$'],['$=(' num2str(winc(1)) ',' num2str(ycent(1)) ')$']};
text(xc+oo(1),yc+oo(2),str,'hori','right', ...
    'fontsize',h.fs(3),'interpreter','latex','verti','bot');


%% LRR 

axes(ax(4));
for j=1:2
    phls(j) = plot(D(j).LRR,'linew',2,'color',colors(j,:)); hold on;
    plot(winc(j),D(j).LRR(winc(j)),'o','linew',2,'color',colors(j,:)); 
end

for j=1:2
    phrr(j) = plot(D(j).relrate,'--','color','k','linew',1); hold on;
    %plot(winc(j),D(j).relrate(winc(j)),'o','linew',2,'color',colors(j,:));
end

%%
tsix = [1 5 4];
axis(ax(tsix),'tight');

set(ax([1 3 4]),'Fontsize',h.fs(end),'XGrid','on','YGrid','on');
set(ax(1),'YTickLabel',[],'XTickLabel',[]);

set(ax(tsix),'Box','off','tickdir','out','fontsize',h.fs(end),h.ticks{:},'XGrid','on','YGrid','on');
axrescale(ax(tsix),0.05, 0.025);


set([H.axd H.axs],'Fontsize',h.fs(end)-4);
set(H.axs(1),'YTickLabel',[]);
set(H.axs(2),'XTickLabel',[]);

ylabel(H.axs(2),[h.mathsym('y') ' index'],'fontsize',h.fs(end));
xlabel(H.axs(1),[h.mathsym('x') ' index'],'fontsize',h.fs(end));

legend(ph,{'$x_i$','$y_j$'},'fontsize',h.fs(3),'location','northeast','interpreter','latex');

xlabel(ax(tsix(end)),'time index','fontsize',h.fs(2));
ylabel(ax(tsix(end-1)), 'rate','fontsize',h.fs(2));
ylabel(ax(tsix(end)), 'LRR','fontsize',h.fs(2));

legend([phls phrr(1)], ...
    {['LRR ' h.mathsym('x') ': \Deltai/\Deltaj'], ...
     ['LRR ' h.mathsym('y') ': \Deltaj/\Deltai'],...
     'rel. process rate'}, ...
    'fontsize',h.fs(2),'interpreter','tex','location','northeast','AutoUpdate','off');

xlabel(ax(3),[h.mathsym('x') ' index (' h.mathsym('i') ')'],'fontsize',h.fs(2));
ylabel(ax(3),[h.mathsym('y') ' index (' h.mathsym('j') ')'],'fontsize',h.fs(2));

labs = {'signals and known process rates',...
    'warping curve',...
    'local relative rate',...
    ['local slope at ' h.mathsym('i') ' = ' num2str(winc(1))]};

axlab = [ax(1) H.axd ax([4 3])];

xoff = [-0.02 -0.02 -0.02 -0.02];

plh = stfig_panlab(axlab,[], ...
    'hori','right','verti','bot','xoff',xoff,'fontsize',h.fs(2),'interp','tex');
stfig_panlab(axlab,labs, ...
    'hori','left','verti','bot','xoff',0.02,'fontsize',h.fs(2),'fontweight','normal');

axes(ax(4));
plot(xlim,[1 1],':','linew',2,'color',[.5 .5 .5]);

%%
h.printfig(mfilename);

end


