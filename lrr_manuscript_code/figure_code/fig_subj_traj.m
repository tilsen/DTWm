function [] = fig_subj_traj()

dbstop if error; close all;
h = tw_helpers;

load([h.datasets_dir 'dataset_CVCarl.mat'],'D');

D.onset = cellfun(@(c)c(1),D.cond);
D.coda = cellfun(@(c)c(end),D.cond);

plot_traj_bycond(D,'onset','coda',h);
h.printfig([mfilename '_byonset']);

plot_traj_bycond(D,'coda','onset',h);
h.printfig([mfilename '_bycoda']);

end

%%
function [] = plot_traj_bycond(D,c0,c1,h)

plot_chans = {'LLy' 'TTy' 'TBy'};
Npch = length(plot_chans);

subjs = unique(D.subj);

info = D.Properties.UserData{end};
chans = info.chans;
t = info.t;

conds = unique(D.(c0));
subconds = unique(D.(c1));
colors = lines(numel(subconds));

condstrs = cellstr(strrep(conds','0',char(216))');
switch(c0)
    case 'onset'        
        subcondstrs = append(condstrs,'.'); 
        condstrs = append('.',condstrs); 

    case 'coda'
        subcondstrs = append('.',condstrs); 
        condstrs = append(condstrs,'.'); 
end

ax = stf([Npch*numel(conds) length(subjs)], ...
    [0.065 0.025 0.01 0.035],[0.01 0.01],'handlearray','matrix','aspect',0.85);

for i=2:Npch
    arrayfun(@(c)shiftposy(c,0.0075*(i-1)),ax(i:Npch:end,:));
end

for s=1:length(subjs)
    for c=1:length(conds)
        for i=1:length(subconds)
            Dx = D(ismember(D.subj,subjs{s}) & D.(c0)==conds(c) & D.(c1)==subconds(i),:);

            for d=1:Npch
                ch_ix = ismember(chans,plot_chans{d});
                axes(ax(Npch*(c-1)+d,s)); %#ok<LAXES> 

                for j=1:height(Dx)
                    ph(s,c,i,d) = plot(t,squeeze(Dx.X(j,ch_ix,:)),'Color',colors(i,:)); hold on;
                end

                if s==1 && i==1
                    ylabel([plot_chans{d}],'fontsize',h.fs(2),'Rotation',0, ...
                        'HorizontalAlignment','right','VerticalAlignment','middle');
                end
            end
        end
        
    end
    title(ax(1,s),subjs{s},'fontsize',h.fs(1));
end

%%
axis(ax,'tight');

legh = legend(squeeze(ph(end,end,:,end)),subcondstrs,'fontsize',h.fs(end),'location','south','NumColumns',3);
legh.Position(2) = ax(Npch,1).Position(2)-legh.Position(4)/2;
legh.Position(1) = ax(Npch,1).Position(1);

set(ax,'xlim',[-0.75 0.75],'XGrid','on','YGrid','on','Box','off','fontsize',h.fs(end)-4);

set(ax,'xtick',[-.5 -.25 0 .25 .5]);
xticklabs = ax(end,1).XTickLabel;
xticklabs([2 4]) = {'',''};
set(ax,'XTickLabel',xticklabs);

set(ax(1:end-1,:),'XTickLabel',[]);
set(ax,'YTickLabel',[],'ytick',-200:5:200);

axrescaley(0.01,ax);

for i=1:Npch
    axi = ax(i:Npch:end,:);
    ylims = arrayfun(@(c){ylim(c)},axi);
    maxrng = max(cellfun(@(c)range(c),ylims(:)));
    arrayfun(@(c)ylim(c,mean(ylim(c))+maxrng*[-1 1]/2),axi);
end


stfig_panlab(ax(1:Npch:end,1),condstrs,'fontsize',h.fs(2), ...
    'yoff',0,'verticalalignment','mid','xoff',-0.15);


end





