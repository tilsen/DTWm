function [] = fig_acoustic_artic_comp()

dbstop if error; close all;
h = tw_helpers;

localrates = {'CVCarl_deltas_100_P05_log','CVCacl_deltas_100_P05_log'};

source = 'empirical';
%source = 'gam';

legstrs = {char(216),'p','t'};
xlab = 'time (s) relative to anchor point';

plot_tlims = [-0.6 0.6];

%%
for i=1:length(localrates)
    switch(source)
        case 'empirical'
            dataset = regexp(localrates{i},'\w+(?=_\w+$)','match','once');
            load([h.lrr_dir 'lrr_' dataset '.mat']);
            conds = unique(G.cond,'stable');
            d.mu = cellfun(@(c){mean(G.lrr_gsd(ismember(G.cond,c),:))},conds);

            all_lrrs = cellfun(@(c){G.lrrs_gsd(ismember(G.cond,c))},conds);
            for j=1:length(conds)
                X = all_lrrs{j};
                X = cell2mat([X{:}]');
                %d.ci{j,1} = norminv(0.975)*std(X)/sqrt(size(X,1));
                d.ci{j,1} = std(X);
            end
            d.cond = conds;
            d.time = repmat(G.t{1}-mean(G.t{1}),length(conds),1);
            D{i} = struct2table(d);
            D{i}.mu = cell2mat(D{i}.mu);
            D{i}.ci = cell2mat(D{i}.ci);

        case 'gam'
            load([h.gamm_dir 'gammpred_' localrates{i} '.mat'],'Ds');
            D{i} = tabindex(Ds,'subj','ALL');
            D{i}.time = D{i}.t-mean(D{i}.t(1,:));
    end
end

ax = stf([3 3],[0.12 0.075 0.01 0.12],[0.01 0.01],'handlearray','matrix');

colors = pastelize([1 0 0; 0 0 1],0.25);
ls = {'-','-'};
lw = [3 3];

conds = unique(D{1}.cond);

for i=1:length(conds)
    [r,c] = ind2sub([3 3],i);
    axes(ax(c,r));
    for j=1:2
        ix = ismember(D{j}.cond,conds(i));
        mu = D{j}.mu(ix,:);
        ci = D{j}.ci(ix,:);
        ci = mu + [-ci; ci];
        t = D{j}.time(ix,:);

        fh(i,j) = fill([t fliplr(t)],[ci(1,:) fliplr(ci(2,:))],colors(j,:),'EdgeColor','none','FaceAlpha',0.25); hold on;
        ph(i,j) = plot(t,mu,'color',colors(j,:),'linew',lw(j),'linestyle',ls{j}); hold on;
    end
end

%%
axis(ax,'tight')
xlim(ax,plot_tlims);
ylim(ax,getlims(ax,'y'));
axrescale(ax,[],0.025);

stfig_panlab(ax(1,:),legstrs,'hori','center','verti','bot','xoff',0.5,'yoff',0,'fontsize',h.fs(1));
stfig_panlab(ax(:,1),legstrs,'hori','center','verti','mid','xoff',-0.2,'yoff',-0.5,'fontsize',h.fs(1));

set(ax,'XGrid','on','YGrid','on','fontsize',h.fs(end),'tickdir','out','box','off','ticklen',0.003*[1 1]);
set(ax(1:end-1,:),'XTickLabel',[]);
set(ax(:,2:end),'YTickLabel',[]);

ylabel(ax(end,1),'LRR-gsd','FontSize',h.fs(2));
xlabel(ax(end,1),xlab,'fontsize',h.fs(3));

legh = legend(ph(1,:),{'articulatory signals','acoustic signals'},'FontSize',h.fs(2),'location','northwest');
shiftposy(legh,0.1);
shiftposx(legh,-0.075);

stbgax;
text(0,0.5,'onset','FontSize',h.fs(1)+4,'Rotation',90,'hori','center','verti','top');
text(ax(1,2).Position*[1 0 0.5 0]',1,'coda','FontSize',h.fs(1)+4,'hori','center','verti','top');


%%
h.printfig(mfilename);

end
