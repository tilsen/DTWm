function [] = fig_anomalous_traj()

dbstop if error; close all;

h = tw_helpers;

%threshold for exclusion of articulatory trajectories
z_thresh = norminv(0.99);

%localrates file:
localrate = 'lrr_CVCarl_deltas_100_P05';

%dataset file:
dataset = 'dataset_CVCarl';

%case for plotting
subj = 'P03';
onset = 'p';
coda = '0';

%channels to plot
plot_chans = {'LLy' 'TTy' 'TBy'};

%%
load([h.datasets_dir dataset '.mat'],'D');
load([h.lrr_dir localrate '.mat'],'G','Gi'); %#ok<*NASGU> 

D.onset = cellfun(@(c)c(1),D.cond);
D.coda = cellfun(@(c)c(end),D.cond);
G.onset = cellfun(@(c)c(1),G.cond);
G.coda = cellfun(@(c)c(end),G.cond);

D = tabindex(D,'subj',subj,'onset',onset,'coda',coda);
G = tabindex(G,'subj',subj,'onset',onset,'coda',coda);

panlabstrs = [strrep([onset 'a' coda],'0',char(216)) ', ' subj];

info = D.Properties.UserData{end};
t = info.t;
chans = info.chans;

ix = 1;
dd = G.dists{ix};
dd(dd==0) = nan;

avgd = nanmean(dd); %#ok<NANMEAN> 
zavgd = zscore(avgd);

avgd_thresh = std(avgd)*z_thresh + mean(avgd);

ix_out = find(zavgd>z_thresh);


%%
ax = stf([1 2 3; 1 2 4],[0.035 0.115 0.01 0.075],[0.065 0.15],'aspect',2.5);

axs = stfig_subaxpos(ax(1),[3 1],[0 0 0 0 .01 .01]);
axs = reshape(axs,[],3);
delete(ax(1));

Nch = length(plot_chans);
gcol = .5*ones(1,3);

for i=1:Nch
    axes(axs(i));
    ch_ix = ismember(chans,plot_chans{i});
    for j=1:height(D)
        ph(i,j) = plot(t,squeeze(D.X(j,ch_ix,:)),'Color',gcol); hold on;
    end

    X = squeeze(D.X(:,ch_ix,:));
    plot(t,mean(X),'color','k','linew',2);
    set(ph(i,ix_out),'color','r','linew',2);
end


%------------------
axes(ax(2));
imh = imagesc(dd);
set(imh,'AlphaData',~isnan(dd));
set(gca,'Color',0.95*[1 1 1],'YDir','normal');

%-------------------
axes(ax(3));
plot(avgd,'.-','color',gcol,'linew',1.5,'markersize',14); hold on;
plot(ix_out,avgd(ix_out),'ro','linew',2);
plot([1 length(avgd)],avgd_thresh*[1 1],'k--','linew',1.5);

%------------------
axes(ax(4));
bh = histogram(avgd,20); hold on;
bh.FaceColor = gcol;
plot(avgd_thresh*[1 1],[0 15],'k--','linew',1.5);
for i=1:length(bh.BinEdges)-1
    if bh.BinEdges(i)>avgd_thresh && bh.BinCounts(i)>0
        fill(bh.BinEdges(i+[0 1 1 0]),[0 0 bh.BinCounts([i i])],'r','FaceAlpha',0.5,'EdgeColor','none');
    end
end

%%

ax = ax(2:end);

axis(axs,'tight');
axrescaley(0.05,axs);

axis(ax(2),'tight');
axrescale(ax(2),0.025,0.05);

set(ax(2),'YGrid','on');
ax(2).XTick(1) = 1;

axis(ax(3),'tight');
axrescale(ax(3),[0 0.05],[0 0.05]);
ax(3).XLim(1) = 0;

set(ax(2:3),'Box','off','YGrid','on','XGrid','on');

set([axs ax],'Fontsize',h.fs(end),'tickdir','out','ticklen',0.003*[1 1],'box','off');

set(axs,'YGrid','on','XGrid','on','XLim',[-.35 .35],'Box','off','xtick',-0.35:.1:.35);
arrayfun(@(c)ylabel(axs(c),[plot_chans{c} ' '], ...
    'fontsize',h.fs(2),'hori','right','rotation',0,'verti','mid'),(1:3));

set(axs(1:end-1),'XTickLabel',[]);
set(axs,'YTickLabel',[]);

stfig_panlab([axs(1) ax], [], 'hori','right','xoff',-0.03,'yoff',-0.01,'verti','bot','fontsize',h.fs(1));

stfig_panlab([axs(1) ax], ...
    {['articulator positions: ' panlabstrs],'pairwise warping distances','average warping distance', 'average distance histogram'}, ...
    'hori','left','xoff',0,'verti','bot','fontsize',h.fs(2),'fontweight','normal');

axs(end).XTickLabelRotation = 0;

xlabel(axs(end),'time (s)','fontsize',h.fs(3));

ax(1).XTick = [1 ax(1).XTick];
ax(1).YTick = [1 ax(1).YTick];
ylabel(ax(1),'trial #','fontsize',h.fs(3));
xlabel(ax(1),'trial #','fontsize',h.fs(3));

ylabel(ax(2),'avg. distance','fontsize',h.fs(3));
xlabel(ax(2),'trial #','fontsize',h.fs(3));

ylabel(ax(3),'count','fontsize',h.fs(3));
xlabel(ax(3),'avg. distance','fontsize',h.fs(3));

%%
h.printfig(mfilename);

end





