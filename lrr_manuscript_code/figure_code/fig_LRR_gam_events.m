function [] = fig_LRR_gam_events()

dbstop if error; close all;
h = tw_helpers;

localrates = 'CVCarl_deltas_100_P05_log';
gamdata = localrates(1:min(regexp(localrates,'_'))-1);

plot_trng = [-0.400 0.600];

load([h.data_dir 'onsetcoda_data.mat']);

ix_keep = ~D.resp_error & ~D.fa_missing & ~isnan(D.V_vel);
D = D(ix_keep,:);

%re-align time coordinate for forced alignment times:
for i=1:height(D)
    P = D.phones{i};
    P.t0 = P.t0-D.V_vel(i);
    P.t1 = P.t1-D.V_vel(i);
    D.phones{i} = P;
end

%%
load([h.gam_dir 'gampred_' localrates '.mat'],'Ds','Dd');

R = load([h.datasets_dir 'summary_' gamdata '.mat']);

Ds.t = Ds.t-mean(Ds.t(1,:));
Dd.t = Dd.t-mean(Dd.t(1,:));

%% collect plotting info

leg_strs = {
    {'/ap/','/pap/'}
    {'/pa/', '/pap/'}};

for i=1:length(leg_strs)
    leg_strs{i} = regexprep(leg_strs{i},'?',char(660));
    leg_strs{i} = regexprep(leg_strs{i},'0',char(216));
end

title_strs = {
    sprintf('onset comparison: %s vs. %s',leg_strs{1}{:}) 
    sprintf('coda comparison: %s vs. %s',leg_strs{2}{:})}; 

title_strs = strrep(title_strs,'{h}','{\fontsize{16}h}');

plot_conds = {...
    {'0_p','p_p'}
    {'p_0','p_p'}};

plot_lms = {...
    {'C1_ons','C1_rel','V_ons','V_trg','C2_ons','C2_trg'},...
    {'C1_ons','C1_rel','V_ons','V_trg','C2_ons','C2_trg'}};

subj = {'ALL','ALL'};

colors = lines(2);

P=[];
for i=1:length(subj)
    P(i).pcs = plot_conds{min(i,end)};
    P(i).lms = plot_lms{i};
    P(i).subj = subj{i};

    P(i).D_lms = tabindex(D,'cond',P(i).pcs);
    P(i).D_data = tabindex(R.D,'cond',P(i).pcs);
    P(i).D_smooths = tabindex(Ds,'cond',P(i).pcs);
    P(i).D_diffs = tabindex(Dd,'cond1',P(i).pcs{1},'cond2',P(i).pcs{2});

    if ~contains(P(i).subj,'ALL')
        P(i).D_lms = tabindex(P(i).D_lms,'subj',P(i).subj);
    end
    
    P(i).D_data = tabindex(P(i).D_data,'subj',P(i).subj);
    P(i).D_smooths = tabindex(P(i).D_smooths,'subj',P(i).subj);
    P(i).D_diffs = tabindex(P(i).D_diffs,'subj',P(i).subj);
    
    P(i).D_smooths.color = colors;
end


%%
axpan = repmat([1 2 3 3 3 4 4 4 5 5 6 6]',1,numel(subj));

for i=2:numel(subj)
    axpan(:,i:end) = axpan(:,i:end) + max(axpan(:))*(i-1);
end

ax = stf(axpan,[0.065 0.065 0.01 0.065],[0.035 0.01],'aspect',1.25);
if size(ax,1)==1, ax=reshape(ax,[],numel(subj)); end
axc = [];
shiftposy(ax(2,:),0.005);

for i=1:length(P)
    hh = plot_events(P(i),ax(:,i));
    set(hh.ax,'xlim',plot_trng);

    HH(i) = hh;
    axc = [axc hh.ax'];

end

%% annotations
bw = 0.005;
xo = -0.225;

axcix = size(axc,1)-[3 2];
for i=1:length(axcix)
    axes(axc(axcix(i))); %#ok<*LAXES> 
    bh = drawbrace([xo -0.5],[xo 1.5],bw,'Color',P(1).D_smooths.color(i,:),'linew',2);
    text(min(bh.XData),mean(bh.YData),leg_strs{1}{i},'hori','right','fontsize',h.fs(3));
end

axcix = size(axc,1)-[1 0];
for i=1:length(axcix)
    axes(axc(axcix(i)));
    bh = drawbrace([xo -0.5],[xo 0.5],bw,'Color',P(1).D_smooths.color(i,:),'linew',2);
    text(min(bh.XData),mean(bh.YData),leg_strs{1}{i},'hori','right','fontsize',h.fs(3));
end

xo = -0.35;
axcix = numel(axc)-[3 2];
for i=1:length(axcix)
    axes(axc(axcix(i))); %#ok<*LAXES> 
    bh = drawbrace([xo -0.5],[xo 1.5],bw,'Color',P(1).D_smooths.color(i,:),'linew',2);
    text(min(bh.XData),mean(bh.YData),leg_strs{2}{i},'hori','right','fontsize',h.fs(3));
end

axcix = numel(axc)-[1 0];
for i=1:length(axcix)
    axes(axc(axcix(i)));
    bh = drawbrace([xo -0.5],[xo 0.5],bw,'Color',P(1).D_smooths.color(i,:),'linew',2);
    text(min(bh.XData),mean(bh.YData),leg_strs{2}{i},'hori','right','fontsize',h.fs(3));
end

axes(axc(end-3,2));
xo = -.175; yo = 0;
phx = plot(xo+[-.03 0],yo*[1 1],'k-'); draw_arrow(phx,1,[0.1 0.01]);
ath(1) = text(phx.XData(1),phx.YData(1),{'constr.' 'onset'},'hori','right');

xo = -.125; yo = 1;
phx = plot(xo+[-.03 0],yo*[1 1],'k-'); draw_arrow(phx,1,[0.1 0.01]);
ath(end+1) = text(phx.XData(1),phx.YData(1),{'vowel' 'onset'},'hori','right');

xo = .025; yo = 0;
phx = plot(xo+[.030 0],yo*[1 1],'k-'); draw_arrow(phx,1,[0.1 0.01]);
ath(end+1) = text(phx.XData(1),phx.YData(1),{'release' 'onset'},'hori','left');

xo = .120; yo = 1;
phx = plot(xo+[.030 0],yo*[1 1],'k-'); draw_arrow(phx,1,[0.1 0.01]);
ath(end+1) = text(phx.XData(1),phx.YData(1),{'target' 'achievement'},'hori','left');

axes(axc(end-1,2));
xo = -.125; yo = 0;
phx = plot(xo+[-.03 0],yo*[1 1],'k-'); draw_arrow(phx,1,[0.1 0.01]);
ath(end+1) = text(phx.XData(1),phx.YData(1),{'segment' 'boundary'},'hori','right');

set(ath,'fontsize',h.fs(end));
text_spacing(ath([1 3]),0.85);
text_spacing(ath([2 4]),0.85);
text_spacing(ath(end),0.95);

%%
axis(ax(1:4,:),'tight');
xlim(ax(1:4,:),plot_trng);

set(axc,'fontsize',h.fs(end));

ylim(ax(1:2,:),getlims(ax(1:2,:),'y'));
ylim(ax(3,:),getlims(ax(3,:),'y'));
ylim(ax(4,:),getlimss(ax(4,:),'y'));
axrescale(ax(1:4,:),0,0.05);

xlabel(axc(end,1:2),'time (relative to anchorpoint)','fontsize',h.fs(3));

ylabel(ax(1,1),'LA','FontSize',h.fs(3),'rotation',0,'hori','right','verti','mid');
ylabel(ax(2,1),'TBy','FontSize',h.fs(3),'rotation',0,'hori','right','verti','mid');
ylabel(ax(3,1),'LRR gsd','fontsize',h.fs(2));
ylabel(ax(4,1),'\Delta LRR gsd','fontsize',h.fs(2));

%
HH(1).legh.String = leg_strs{1};
HH(2).legh.String = leg_strs{2};
set([HH.legh],'Fontsize',h.fs(3));

%
labels = {'(a)' '(b)' '(c)' '(d)' '(e)' '\alpha' '\beta'};
axi = [1 1 1 1 1 2 2];
axj = [1 1 1 1 1 1 2];
times = [-0.200   0    0.125    0.250  0.35 -0.155 0.165];
yposn = [1 .5 .6 .5 .65 0.425 0.375];
for i=1:length(labels)
    axes(axc(axi(i)+2,axj(i)));
    lth(i) = text(times(i),min(ylim)+diff(ylim)*yposn(i),labels{i}, ...
        'fontsize',h.fs(3),'hori','center','verti','mid');
end
set(lth(end-1:end),'fontsize',h.fs(1));

%anchorpoint lines:
arrayfun(@(c)set(c,'XData',[0 0],'Ydata',ylim(c.Parent),'Linestyle',':','color',0.5*ones(1,3),'linew',2), ...
    [HH.anchline]);

stfig_panlab(ax(1,:),[],'xoff',0,'hori','right');
stfig_panlab(ax(1,:),title_strs,'hori','left','xoff',0.025,'fontweight','normal');

stbgax;
text(0.065,axc(end-3,1).Position(2),{'gestural','intervals'},...
    'rotation',90,'verti','bot','hori','center','fontsize',h.fs(2));

text(0.065,axc(end-1,1).Position(2),{'phones'},...
    'rotation',90,'verti','bot','hori','center','fontsize',h.fs(2));

set(axc(:,2),'YTickLabel',[]);
axc(end,2).XTickLabel{1} = '';

stfig_panlab(axc(4,1),strcat(leg_strs{1}{1},' -',[' ' leg_strs{1}{2}]), 'hori','left','xoff',0.01, ...
    'yoff',-1,'verti','bot','fontw','normal');
stfig_panlab(axc(4,2),strcat(leg_strs{2}{1},' -',[' ' leg_strs{2}{2}]), 'hori','left','xoff',0.01, ...
    'yoff',-1,'verti','bot','fontw','normal');

set(axc(end-3:end,:),'GridAlpha',1,'GridColor',[0 0 0]);

F = tw_extract_lrrgsd_features(localrates);
feats = {'alpha_vel','beta_vel'};

for i=1:numel(subj)
    axes(axc(4,i));
    Fx = tabindex(F,'subj',subj{i},'cond1',plot_conds{i}{1},'cond2',plot_conds{i}{2});
    featval = Fx.(feats{i});
    plot(featval*[1 1],[0 max(ylim)],'r-','linew',2);
end

%%

h.printfig(mfilename);

end



%%
function [h] = plot_events(PP,ax)

sefac = norminv(0.95);

Dx = PP.D_lms;
Ds = PP.D_smooths;
Dd = PP.D_diffs;
Dr = PP.D_data;
conds = PP.pcs;
lms = PP.lms;

dens_pnts = -0.75:0.001:0.750;

for i=1:length(conds)
    dx = tabindex(Dx,'cond',conds{i});
    X = table2array(dx(:,lms))-dx.V_vel;
    ix_valid = ~all(isnan(X));
    P(i).lmx = lms(ix_valid);
    P(i).X = X(:,ix_valid);    
    for j=1:length(P(i).lmx)   
        [P(i).dens(j,:),x] = ksdensity(P(i).X(:,j),dens_pnts);
        [P(i).fmax(j),P(i).fmax_ix(j)] = max(P(i).dens(j,:));
        P(i).fmax_t(j) = x(P(i).fmax_ix(j));
        %P(i).fmax_t(j) = nanmean(P(i).X(:,j)); %#ok<NANMEAN>      
    end
    P(i).C_ix = contains(P(i).lmx,'C');
    P(i).V_ix = contains(P(i).lmx,'V');
end

%articulator trajectories
info = Dr.Properties.UserData{1};
chix = @(ch)ismember(info.kin_chans,ch);
for i=1:length(conds)

    axes(ax(1)); h.anchline(1) = plot(nan,nan); hold on;
    axes(ax(2)); h.anchline(2) = plot(nan,nan); hold on;

    rx = tabindex(Dr,'cond',conds{i});
    T=[];
    if contains(rx.cond{1},'p')
        LL = squeeze(rx.mu(1,chix({'LLx' 'LLy'}),:));
        UL = squeeze(rx.mu(1,chix({'ULx' 'ULy'}),:));
        x = -sqrt(sum((LL-UL).^2));
        T(end+1).lab = 'LA';
        T(end).x = x-mean(x);
        T(end).panix = 1;
    end
    if contains(rx.cond{1},'t')
        x = squeeze(rx.mu(1,chix({'TTy'}),:))';
        T(end+1).lab = 'TTy';
        T(end).x = x-mean(x);        
        T(end).panix = 1;
    end
    x = squeeze(rx.mu(1,chix({'TBy'}),:))';
    T(end+1).lab = 'TBy';
    T(end).x = x-mean(x);
    T(end).panix = 2;
    
    T = struct2table(T);

    for j=1:height(T)
        h.trph(i,j) = plot(rx.t(1,:),T.x(j,:), ...
            '-','color',Ds.color(i,:),'linew',2,'parent',ax(T.panix(j))); hold on;
        h.trchans{i,j} = T.lab{j};
    end

end

%LRR-gsd
axes(ax(3));
for i=1:height(Ds)
    h.anchline(3) = plot(nan,nan); hold on;
    [fh(i),ph(i)] = traj_plot(Ds.mu(i,:),sefac*Ds.ci(i,:),Ds.t(i,:),'Color',Ds.color(i,:)); hold(gca,'on');
end
h.legh = legend(ph,conds,'location','southwest');

%LRR-gsd diff
axes(ax(4));
h.anchline(4) = plot(nan,nan); hold on;
plot(minmax(dens_pnts),[0 0],'k--','linew',2); hold on;
[dfh,dph] = traj_plot(Dd.mu(1,:),sefac*Dd.ci(1,:),Dd.t(1,:));


%kinematic landmarks
axes(ax(5));
axk = stfig_subaxpos(gca,[2 1],[0 0 0 0 0 0.0025]);
for i=1:length(conds)

    axes(axk(i));
    h.anchline(end+1) = plot([0 0],[nan nan]); hold on;

    color = Ds.color(i,:);

    pp = P(i);
    f = zeros(2,length(pp.dens));
    f(1,:) = sum(pp.dens(pp.C_ix,:),1);
    f(2,:) = sum(pp.dens(pp.V_ix,:),1);

    kinfh(i) = imagesc(dens_pnts,[0 1],f,'AlphaData',0.95*ones(size(f))); hold on;
    cmap = ST_colormap(1000,[1 1 1],color);
    colormap(gca,cmap);

    for j=1:length(pp.lmx)
        oy = double(pp.V_ix(j));
        yy = 0.5*[-1 1]+oy;
        tstr = regexp(pp.lmx{j},'^\w+(?=_)','match','once');
        switch(contains(pp.lmx{j},{'trg' 'rel'}))
            case 1
                plot(pp.fmax_t(j)*[1 1],yy,':','color',color,'linew',2); hold on;
                %text(pp.fmax_t(j),oy,tstr,'fontsize',16,'hori','right','interp','none');
            otherwise
                plot(pp.fmax_t(j)*[1 1],yy,'-','color',color,'linew',2); hold on;
                text(pp.fmax_t(j),oy,[' ' tstr],'fontsize',16,'hori','left','interp','none');
        end

        switch(contains(pp.lmx{j},'ons'))
            case 1
                switch(contains(pp.lmx{j},'C1'))
                    case 1
                        ix_trg = find(contains(pp.lmx,strrep(pp.lmx{j},'ons','rel')));
                        xx = pp.fmax_t([j ix_trg]);
                        gacth(i,j) = plot(xx([1 2 2 1 1]),yy([1 1 2 2 1]),'-','color',color,'linew',2);                        
                    otherwise
                        ix_trg = find(contains(pp.lmx,strrep(pp.lmx{j},'ons','trg')));
                        xx = pp.fmax_t([j ix_trg]);
                        gacth(i,j) = plot(xx([2 1 1 2]),yy([1 1 2 2]),'-','color',color,'linew',2);
                        gacthe(i,j) = plot(xx([2 2]),yy([1 2]),':','color',color,'linew',2);                        
                end

        end
    end


end

f=[];

%segmentation
axes(ax(6));
colors = lines(2);
axs = stfig_subaxpos(gca,[2 1],[0 0 0 0 0 0.0025]);
for i=1:length(conds)

    axes(axs(i));
    h.anchline(end+1) = plot(nan,nan); hold on;

    dx = tabindex(Dx,'cond',conds{i});
    P = vertcat(dx.phones{:});
    P = tabindex(P,'~lab','SIL');
    segs = unique(P.lab,'stable');
    %colors = hsv(numel(segs));

    for j=1:length(segs)
        X = P.t1(ismember(P.lab,segs{j}));
        [f(j,:),x] = ksdensity(X,dens_pnts);
    end

    ff = sum(f);
    segfh(i) = imagesc(dens_pnts,0,ff,'AlphaData',0.95*ones(size(ff))); hold on;
    cmap = ST_colormap(1000,[1 1 1],colors(i,:));
    colormap(gca,cmap);

    for j=1:length(segs)
        [fmax(j),fmax_ix(j)] = max(f(j,:));
    end
    tsegs = x(fmax_ix);
    for j=1:length(segs)-1
        plot(tsegs(j+[0 1 1 0 0]),0.5*[-1 -1 1 1 -1],'-','color',colors(i,:),'linew',2); hold on;       
        text(mean(tsegs(j+(0:1))),0,strrep(segs{j+1},'q',['{\fontname{Cambria}' char(660) '}']), ...
            'fontsize',18,'hori','center');
    end


end


%%
delete(ax([5:6]));

ax = [ax(1:4)' axk axs];

set(axk,'YDir','reverse');
axis([axk axs],'tight');
set([axk axs],'box','off','tickdir','out','ycolor','none');
set(ax([5 7]),'xcolor','none');

set(ax(1:4),'YGrid','on');
set(ax(1:2),'YTick',[]);

set(ax,'box','off','tickdir','out','ticklen',0.003*[1 1],'xgrid','on');

set(ax(1:end-1),'XTickLabel',[]);

h.ax = ax;


end


