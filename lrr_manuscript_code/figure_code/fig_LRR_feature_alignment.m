function [] = fig_LRR_feature_alignment()

dbstop if error; close all;
h = tw_helpers;

localrates = 'CVCarl_deltas_100_P05_log';

F = tw_extract_lrrgsd_features(localrates);
load([h.data_dir 'onsetcoda_data.mat']);

%exclude P04-alpha feature values
F.alpha_vel(ismember(F.subj,'P04') & F.ax==true) = nan;

%exclude bad data
ix_keep = ~D.resp_error & ~D.fa_missing & ~isnan(D.V_vel);
D = D(ix_keep,:);

%re-align time coordinates:
for i=1:height(D)
    P = D.phones{i};
    P.t0 = P.t0-D.V_vel(i);
    P.t1 = P.t1-D.V_vel(i);
    D.phones{i} = P;
end

lms = {'C1_ons','C1_rel','V_ons','V_trg','C2_ons','C2_trg'};
for i=1:length(lms)
    D.(lms{i}) = D.(lms{i}) - D.V_vel;
end

%generalize segmentation
D.C1_t0 = cellfun(@(c)c.t0(ismember(c.lab,{'q' 'p1' 't1'})),D.phones);
D.C1_t1 = cellfun(@(c)c.t1(ismember(c.lab,{'q' 'p1' 't1'})),D.phones);
D.V_t0 = cellfun(@(c)c.t0(ismember(c.lab,{'aa'})),D.phones);
D.V_t1 = cellfun(@(c)c.t1(ismember(c.lab,{'aa'})),D.phones);

D.C2_t0 = nan(height(D),1);
D.C2_t1 = nan(height(D),1);
D.C2_t0(ismember(D.coda,{'p' 't'})) = cellfun(@(c)c.t0(ismember(c.lab,{'p2' 't2'})),D.phones(ismember(D.coda,{'p' 't'})));
D.C2_t1(ismember(D.coda,{'p' 't'})) = cellfun(@(c)c.t1(ismember(c.lab,{'p2' 't2'})),D.phones(ismember(D.coda,{'p' 't'})));

%%

lms = [lms {'C1_t0' 'V_t0' 'C2_t0' 'C2_t1'}];

G = grpstats(D,{'subj' 'cond'},{@(x)nanmean(x)},'DataVars',lms); %#ok<NANMEAN> 
G.Properties.VariableNames = regexprep(G.Properties.VariableNames,'Fun1_','');

subjs = unique(D.subj);

comps_onset = {
    '0_p' 'p_p';
    '0_t' 'p_t';
    '0_0' 'p_0';
    '0_p' 't_p';
    '0_t' 't_t';
    '0_0' 't_0'};

comps_coda = cellfun(@(c){fliplr(c)},comps_onset);

C = {comps_onset,comps_coda};
feats = {'alpha_vel' 'beta_vel'};
E = {{'C1_ons' 'C1_t0'} {'C2_ons' 'C2_t0'}};


%%
ax = stf([numel(subjs)+1 numel(C)],[0.035 0.075 0.01 0.05],[0.05 0.005],'handlearray','matrix');

vlh = [];
ecol = pastelize([0 0 1; 1 0 0],0.5);

for i=1:length(C)
    
    comps = C{i};
    events = E{i};
    
    for k=1:size(comps,1)
        
        axes(ax(k,i));
        vlh(end+1) = plot(nan,nan,'k--'); hold on;

        for j=1:length(subjs)

            [~,ix] = tabindex(F,'subj',subjs{j},'cond1',comps(k,1),'cond2',comps(k,2));
            featx = F.(feats{i})(ix);
            
            Gx = [tabindex(G,'subj',subjs{j},'cond',comps(k,1)); tabindex(G,'subj',subjs{j},'cond',comps(k,2))];
            y = table2array(Gx(end,events))-featx;
            
            ms = [8 10];
            for m=1:length(events)
                switch(contains(events{m},{'t0' 't1'})) %segmentation
                    case 1
                        ls = 's';
                        color = ecol(m,:);
                    otherwise %gestural landmark
                        ls = 'o';
                        color = ecol(m,:);
                end

                phe(i,k,j,m) = plot(y(m),1-j,ls, ...
                    'color','k','markerfacecolor',color, ...
                    'markersize',ms(m),'linew',1); hold on;
            end

            Y(i,k,j,:) = y;
        end
    end
end

%%
Y = permute(Y,[1 4 2 3]);
Y = reshape(Y,size(Y,1),size(Y,2),[]);

for i=1:length(C)
    axes(ax(end,i));

    vlh(end+1) = plot(nan,nan,'k--'); hold on;

    y = squeeze(Y(i,:,:));

    for j=1:size(y,1)
        [dens,pnts] = ksdensity(y(j,:),'numpoints',501);
        [hd(i,j),ph(i,j)] = dens_plot(dens,pnts,'color',ecol(j,:),'outline',true); hold on;
    end

end

%%

axis(ax,'tight');
set(ax(:,1),'xlim',getlimss(ax(:,1),'x'));
set(ax(:,2),'xlim',getlimss(ax(:,2),'x'));
axrescale(ax(1:end-1,:),0.10,0.10);
axrescale(ax(end,:),0.10,[0 0.05]);

arrayfun(@(c)set(c,'YData',ylim(get(c,'parent')),'XData',[0 0],'linewidth',2),vlh);

set(ax,'Box','off','yticklabel',[],'tickdir','out','ticklen',0.003*[1 1], ...
    'fontsize',h.fs(end),'XGrid','on');

cons = regexprep(comps_onset,'^(\w{1})','{\\bf{$1}}');
ccod = regexprep(comps_coda,'(\w{1})$','{\\bf{$1}}');

comps = [cons; ccod];

%comps = cellfun(@(c){['/' c '/']},comps);
comps = regexprep(comps,'0',char(216));
comps = regexprep(comps,'_','a');

panlabs = cellfun(@(c,d){[c ' - ' d]},comps(:,1),comps(:,2));

stfig_panlab(ax(1:end-1,:),panlabs, ...
    'hori','left','xoff',0.01,'yoff',0,'verti','top','fontsize',h.fs(1),'fontweight','normal');

set(ax(1:end-1,:),'XTickLabel',[],'Ytick',-5:0);
xlabel(ax(end,1),'time (s) relative to \alpha','fontsize',h.fs(3));
xlabel(ax(end,2),'time (s) relative to \beta','fontsize',h.fs(3));
ylabel(ax(end-1,1),'participants','fontsize',h.fs(end));

stfig_panlab(ax(1,:),[],'hori','right','xoff',0,'fontsize',h.fs(1));
stfig_panlab(ax(1,:),{'event alignment to \alpha','event alignment to \beta'}, ...
    'hori','left','xoff',0.025,'fontweight','normal','fontsize',h.fs(1));

set(ax(end,:),'YColor','none');

objs = squeeze(phe(1,6,1,:));
legh(1) = legend(objs,{'gest-init. (C1)' 'bound. (C1)'},'location','southeast','fontsize',h.fs(2));
shiftposx(legh,0.025);

objs = squeeze(phe(2,6,1,:));
legh(2) = legend(objs,{'gest-init. (C2)' 'bound. (C2)'},'location','southeast','fontsize',h.fs(2));

%%
h.printfig(mfilename);

end







