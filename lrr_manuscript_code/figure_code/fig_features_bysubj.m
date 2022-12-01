function [] = fig_features_bysubj()

dbstop if error; close all;
h = tw_helpers;

setnames = {'onsets','codas'};

for i=1:length(setnames)
    plot_features(setnames{i},h);
    close all;
end

end

%%
function [] = plot_features(setname,h)

localrates = 'CVCarl_deltas_100_P05_log';

load([h.gam_dir 'gampred_' localrates '.mat'],'Ds','Dd');

F = tw_extract_lrrgsd_features(localrates);

Ds.t = Ds.t-mean(Ds.t(1,:));
Dd.t = Dd.t-mean(Dd.t(1,:));

plot_tlims = [-0.4 0.6];

%%
subjs = unique(Dd.subj);

plot_conds = {
    '0_p' 'p_p';
    '0_t' 'p_t';
    '0_0' 'p_0';
    '0_p' 't_p';
    '0_t' 't_t';
    '0_0' 't_0'};

switch(setname)
    case 'codas'
        plot_conds = cellfun(@(c){fliplr(c)},plot_conds);
end


ax = stf([numel(subjs) size(plot_conds,1)],[0.035 0.055 0.01 0.035],[0.005 0.0025], ...
    'handlearray','matrix');

for a=1:size(ax,1)
    for b=1:size(ax,2)
        axc(a,b,:) = stfig_subaxpos(ax(a,b),[2 1],[0 0 0 0 0 0]);
        delete(ax(a,b));
    end
end

cifcn = @(dx)[dx.mu(1,:)-dx.ci(1,:); dx.mu(1,:)+dx.ci(1,:)];
cifcn_se = @(dx,fac)[dx.mu(1,:)-fac*dx.ci(1,:); dx.mu(1,:)+fac*dx.ci(1,:)];

colors = lines(2);

h.phf = [];
feats = F.Properties.VariableNames(contains(F.Properties.VariableNames,{'alpha' 'beta'}));

for i=1:size(plot_conds,1)
    for j=1:length(subjs)

        axes(axc(j,i,1));
        Dsj = tabindex(Ds,'subj',subjs{j},'cond',plot_conds(i,:));

        for k=1:2
            ci = cifcn_se(Dsj(k,:),norminv(0.975));            
            traj_plot(Dsj.mu(k,:),ci,Dsj.t(k,:),'color',colors(k,:)); hold on;
        end

        axes(axc(j,i,2));
        Ddj = tabindex(Dd,'subj',subjs{j},'cond1',plot_conds(i,1),'cond2',plot_conds(i,2));
        ci = cifcn_se(Ddj,1);

        plot(minmax(Ddj.t(1,:)),[0 0],'k--'); hold on;
        traj_plot(Ddj.mu(1,:),ci,Ddj.t(1,:),'color',[0 0 0]);

        Fj = tabindex(F,'subj',subjs{j},'cond1',plot_conds(i,1),'cond2',plot_conds(i,2));
        
        if ~isempty(Fj)
            x = table2array(Fj(:,feats));
            for m=1:length(x)
                if isnan(x(m)), continue; end
                yix = find(Ddj.t(1,:)>=x(m),1,'first');
                h.phf(end+1) = plot(x(m)*[1 1],ci(:,yix)','r-','linew',2);
            end
        end
    end
end

axis(axc,'tight');

set(axc(:,:,1),'ylim',getlims(axc(:,:,1),'y'));
set(axc(:,:,2),'ylim',getlims(axc(:,:,2),'y'));

axrescale(axc,[],0.025);
set(axc,'YGrid','on','XGrid','on','Ticklen',0.003*[1 1],'XLim',plot_tlims,'fontsize',h.fs(end));

set(axc(1:end-1,:,:),'XTickLabel',[]);
set(axc(1:end,:,1),'XTickLabel',[]);
set(axc,'YTickLabel',[]);

xlabel(axc(end,1,2),'time (relative to anchorpoint)','FontSize',h.fs(end));

stfig_panlab(axc(:,1,1),subjs,'xoff',-0.01,'yoff',-0.5,'hori','right','fontsize',h.fs(3));
stfig_panlab(axc(1,:,1),cellfun(@(c,d){[c ' - ' d]},plot_conds(:,1),plot_conds(:,2)), ...
    'xoff',0.5,'yoff',0,'hori','center','fontsize',h.fs(3),'verti','bot','interpreter','none');

arrayfun(@(c)set(c,'YData',ylim(get(c,'Parent'))),h.phf);

%%
h.printfig(strjoin({mfilename setname},'_'));

end

