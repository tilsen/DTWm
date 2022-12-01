function [] = fig_LRR_CVC_both()

dbstop if error; close all;
h = tw_helpers;

localrates = 'lrr_CVCarl_deltas_100_P05';
%localrates = 'lrr_CVCacl_deltas_100_P05';

%%
load([h.lrr_dir localrates '.mat'],'G','Gi');

subjs = unique(G.subj);
Nc = length(unique(G.onset));
Ns = length(subjs);


%%
ax = stf([1+Ns 2*Nc], [0.135 0.065 0.01 0.125],[0.005 0.005], ...
    'handlearray','matrix');

shiftposx(ax(:,1:end/2),-0.065);
adjheight(ax(1,:),0.025);

axm = ax(1,:);
axs = ax(2:end,:);


groupf = reshape(repmat({'codaix' 'onsetix'},3,1),1,[]);
plotf = reshape(repmat({'onsetix' 'codaix'},3,1),1,[]);
segixs = repmat(1:Nc,1,2);
seglabs = {char(216),'p','t'};

colors = lines(Nc);

for i=1:size(ax,2)
    
    axes(axm(i));
    Gx = tabindex(G,groupf{i},segixs(i));

    for j=1:Nc
        Gy = tabindex(Gx,plotf{i},j);

        x = mean(Gy.lrr_gsd);
        t = Gy.t{1};
        ph(i,j) = plot(t,x,'color',colors(j,:),'linew',3); hold on;

    end

    for j=1:Ns
        axes(axs(j,i));
        for k=1:Nc
            Gy = tabindex(Gx,'subjix',j,plotf{i},k);   
            x = Gy.lrr_gsd;
            t = Gy.t{1};
            phs(k,i,j) = plot(t,x,'color',colors(k,:),'linew',2); hold on;
        end
    end

end

%%
axis(ax,'tight');
set(ax,'ylim',getlims(ax,'y'));
axrescale(ax,0.01,0.025);

set(ax,'XGrid','on','YGrid','on','fontsize',h.fs(end));

set(ax(:,[2 3 5 6]),'YTickLabel',[]);
set(ax(1:end-1,:),'XTickLabel',[]);

sl = seglabs(segixs);
sl(1:3) = append('coda ',sl(1:3));
sl(4:end) = append('onset ',sl(4:end));

stfig_panlab(ax(1,:),sl,'hori','center','verti','bot', ...
    'fontsize',h.fs(2),'xoff',0.5,'yoff',0,'fontweight','normal');

stfig_panlab(ax(:,1),[{'AVG.'}; subjs],'hori','left','verti','mid', ...
    'fontsize',h.fs(2),'xoff',-0.45,'yoff',-0.5,'fontweight','normal');

legh(1) = legend(ph(3,:),seglabs,'location','north','fontsize',h.fs(2),'NumColumns',Nc);
legh(2) = legend(ph(6,:),seglabs,'location','north','fontsize',h.fs(2),'NumColumns',Nc);

stbgax;
text(0.01,0.975,'A) LRR-gsd comparison by onset','FontSize',h.fs(2)+2,'fontweight','bold');
text(0.55,0.975,'B) LRR-gsd comparison by coda','FontSize',h.fs(2)+2,'fontweight','bold');

text(axs(end,2).Position*[1 0 0.5 0]',0.001,'time (s) relative to anchorpoint', ...
    'verti','bot','hori','center','fontsize',h.fs(3));

text(axs(end,5).Position*[1 0 0.5 0]',0.001,'time (s) relative to anchorpoint', ...
    'verti','bot','hori','center','fontsize',h.fs(3));

strs = {'onsets: ','codas: '};

shiftposy(legh,0.095);
for i=1:2
    text(legh(i).Position(1),legh(i).Position*[0 1 0 0.5]',strs{i},'hori','right','fontsize',h.fs(2));
end

ylabel(ax(1,4),'LRR-gsd','fontsize',h.fs(2));

PP = [-0.6 -0.30; 
    -0.25 0.30; 
    0.35 0.6];

yo = 2;

phaselabs = {'(i)','(ii)','(iii)'};
axes(ax(1,1));
for i=1:size(PP,1)
    drawbrace([PP(i,1) yo],[PP(i,2) yo],0.005,'color','k','linew',1.5);
    text(mean(PP(i,:)),yo+0.1,phaselabs{i},'hori','center','verti','bot','fontsize',h.fs(3),'fontweight','bold');
end


%%

h.printfig([mfilename '_' strrep(localrates,'localrates_','')]);

end

