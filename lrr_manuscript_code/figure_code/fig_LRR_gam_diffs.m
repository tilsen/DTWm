function [] = fig_LRR_gam_diffs()

dbstop if error; close all;
h = tw_helpers;

localrates = 'CVCarl_deltas_100_P05_log';

plot_trng = [-0.650 0.650];

legstrs = {char(216),'p','t'};
segments = {'0','p','t'};
xlab = 'time (s) relative to anchor point';

sefac = norminv(0.975);

%%
load([h.gam_dir 'gampred_' localrates '.mat'],'Ds','Dd');
Ds.time = Ds.t-mean(Ds.t,2);
Dd.time = Dd.t-mean(Dd.t,2);

conds = unique(Ds.cond,'stable');
Ds = tabindex(Ds,'subj','ALL');
Dd = tabindex(Dd,'subj','ALL');
Ds.cond = grp2idx(Ds.cond);

onsets = cellfun(@(c){c(1)},conds);
codas = cellfun(@(c){c(end)},conds);

ax = stf([4 6],[0.085 0.065 0.01 0.155],[0.005 0.0035],'handlearray','matrix');

shiftposx(ax(:,1:3),-0.025);
for i=1:size(ax,2)
    ax(1,i).Position(4) = ax(1,i).Position(4)+0.075;
end

colors = repmat(lines(3),3,1);

ls = {'-',':'};
lw = [3 3];

for c=1:2 %onset, coda comparisons
    for a=1:length(segments) %fixed onset or coda

        axes(ax(1,a + (c-1)*3));

        switch(c)
            case 1
                for i=1:length(segments)
                    cix(i) = find(ismember(onsets,segments{i}) & ismember(codas,segments{a}));
                end

            case 2
                for i=1:length(segments)
                    cix(i) = find(ismember(codas,segments{i}) & ismember(onsets,segments{a}));
                end
        end

        for i=1:length(cix)
            y = Ds.mu(Ds.cond==cix(i),:);
            ci = Ds.ci(Ds.cond==cix(i),:);
            ci = y + sefac*[-ci; ci];
            t = Ds.time(Ds.cond==cix(i),:);
            color = colors(i,:);
            fh(c,a,i) = fill([t fliplr(t)],[ci(1,:) fliplr(ci(2,:))],color,'EdgeColor','none','FaceAlpha',0.25); hold on;
            ph(c,a,i) = plot(t,y,'color',color,'linew',lw(1),'linestyle',ls{1}); hold on;

            legstrs{c,a,i} = ['.' strrep(conds{cix(i)},'0',char(216)) '.'];
        end

        legend(squeeze(ph(c,a,:)),squeeze(legstrs(c,a,:)), ...
            'Interpreter','none','fontsize',h.fs(end),'location','southwest','NumColumns',1,'autoupdate',false);

        CC = [1 2; 1 3; 2 3];

        for i=1:size(CC,1)

            rowix = 1+i;
            colix = a + (c-1)*3;
            axes(ax(rowix,colix));

            switch(c)
                case 1 %fixed coda, onset varies
                    cix1 = find(ismember(onsets,segments{CC(i,1)}) & ismember(codas,segments{a}));
                    cix2 = find(ismember(onsets,segments{CC(i,2)}) & ismember(codas,segments{a}));
                    diffstr{rowix,colix} = ['.' segments{CC(i,1)} ' - ' '.' segments{CC(i,2)}];

                case 2 %fixed onset, coda varies
                    cix1 = find(ismember(codas,segments{CC(i,1)}) & ismember(onsets,segments{a}));
                    cix2 = find(ismember(codas,segments{CC(i,2)}) & ismember(onsets,segments{a}));
                    diffstr{rowix,colix} = [segments{CC(i,1)} '.' ' - ' segments{CC(i,2)} '.'];

            end

            [~,ixd] = tabindex(Dd,'cond',sprintf('%i_%i',cix1,cix2));
            y = Dd.mu(ixd,:);
            ci = Dd.ci(ixd,:);
            ci = y + sefac*[-ci; ci];
            t = Dd.time(ixd,:);
            color = [0 0 0];

            plot(minmax(t),[0 0],'k:','linew',2); hold on;

            fill([t fliplr(t)],[ci(1,:) fliplr(ci(2,:))],color,'EdgeColor','none','FaceAlpha',0.25); hold on;
            plot(t,y,'color',color,'linew',lw(1),'linestyle',ls{1}); hold on;

            [~,maxdiffix] = find(abs(y)==max(abs(y)));
            max_diff_t(c,a,i) = t(maxdiffix);

            set(gca,'UserData',t(maxdiffix));

        end
    end
end

%%
axis(ax,'tight')

axs = ax(1,:);
axd = ax(2:end,:);

ylim(axs,getlims(axs,'y'));
ylim(axd,getlimss(axd,'y'));

%getlims(axs,'y')
set(axs,'ytick',0:0.2:10);
set(axd,'ytick',-0.4:0.2:.4);
set(ax,'XTick',-.5:0.25:.5);
set(ax,'XLim',plot_trng);

axrescale(ax,0.01,0.025);

set(ax,'YGrid','on','XGrid','on','ticklen',0.003*[1 1],'fontsize',h.fs(end)-2,'box','off');

set(ax(1:end-1,:),'XTickLabel',[]);
set(ax(:,2:end),'YTickLabel',[]);

axpanlabs = regexprep(segments,'0',char(216));
plprops = {'hori','center','verti','bot','xoff',0.5,'yoff',0,'fontsize',h.fs(3),'fontweight','normal'};
plprops1 = {'hori','left','verti','bot','xoff',0.01,'yoff',0,'fontsize',h.fs(3),'fontweight','normal'};
stfig_panlab(ax(1,1:3),cellfun(@(c){['coda ' c]},axpanlabs),plprops1{:});
stfig_panlab(ax(1,4:end),cellfun(@(c){['onset ' c]},axpanlabs),plprops1{:});

diffstr = diffstr(2:end,:);
diffstr = strrep(diffstr,'0',char(216));
plh = stfig_panlab(axd(:),diffstr(:),'verti','top','hori','left','xoff',0.01,'yoff',0,'fontsize',h.fs(3));
%set(plh,'BackGroundColor','w');

ylabel(ax(1,1),'LRR-gsd','FontSize',h.fs(3));
ylabel(ax(2:end,1),'\Delta LRR-gsd','FontSize',h.fs(3));
xlabel(ax(end,[2 5]),xlab,'fontsize',h.fs(3));

%arrayfun(@(c)plot([0 0],ylim(c),'k:','linew',1.5,'parent',c),ax);

stbgax;
text(ax(1,1).Position(1)-0.02,0.98,'A) comparisons by onset','fontsize',h.fs(2),'fontweight','bold');
text(ax(1,4).Position(1)-0.02,0.98,'B) comparisons by coda','fontsize',h.fs(2),'FontWeight','bold');

yo = [-0.25 -0.01];
for j=1:6
    for k=1:2
        axes(axd(k,j));
        xo = axd(k,j).UserData;
        arrow([xo yo(1)],[xo yo(2)],'length',8,'tipangle',30,'linew',2);
    end
end

%%
h.printfig(mfilename);

end
