function [] = fig_edge_effect()

dbstop if error; close all;
h = tw_helpers;

step_patterns = {'symmetric2','symmetricP05'};

S{1} = test_signals('edge_effect_example');

D(1).X{1} = S{1}.X{1};
D(1).X{2} = S{1}.X{2};
D = repmat(D,2,1);

%
for i=1:length(D)
    [D(i).map,D(i).distance,D(i).map_info] = dtwm(D(i).X{1},D(i).X{2},'step_pattern',step_patterns{i});

    [D(i).Xa(1,:),D(i).Xa(2,:)] = align_signals(D(i).map,D(i).X{1},D(i).X{2},'aligntype','expansive');
end

%%
axpan = [1 3; 2 3; 4 6; 5 6];
axpan = axpan(:,[1 1 1 2 2]);
ax = stf(axpan,[0.025 0.065 0.035 0.05],[0.075 0.075],'aspect',1.5);
ax = reshape(ax,[],2)';

for j=1:2
    shiftposy(ax(j,[1 2]),[-0.035 0.035]);
end

colors = lines(2);
lw = [2 3];
ls = {'-',':'};

yo = 0.0;

for i=1:length(D)

    axes(ax(i,1)); %#ok<*LAXES> 
    for j=1:2
        ph(i,j) = plot(D(i).X{j}+yo*j,ls{j},'color',colors(j,:),'linew',lw(j)); hold on;
    end

    if i==1
        ixs = [find(D(i).X{2}>=D(i).X{1}(end),1,'first') length(D(i).X{1})];
        plot(ixs,[D(i).X{2}(ixs(1)) D(i).X{1}(end)],'ko-','linew',2);
    end

    axes(ax(i,2));
    for j=1:2
        pha(i,j) = plot(D(i).Xa(j,:)+yo*j,ls{j},'color',colors(j,:),'linew',lw(j)); hold on;
    end    

    axes(ax(i,3));
    hd{i} = plot_distance(D(i).X{1},D(i).X{2},D(i).map);

    set(hd{i}.mh,'LineWidth',3,'Color',pastelize([1 0 0],0.5));

    xlim(xlim+[-1 1]*3);
    ylim(ylim+[-1 1]*3);
       
end

%%

axs = ax(:,[1 2]);
axis(axs,'tight');
axrescaley(axs,0.05);
set(axs,'XLim',[1 length(D(1).Xa)+1]);
set(axs,'XTick',[1 200:200:length(D(1).Xa)+1]);
axrescalex(axs,0.01);

set(ax(:,1),'XTickLabel',[]);
set(axs,'YTickLabel',[],'Box','off');

set(ax,'TickDir','out','ticklen',0.003*[1 1],'Fontsize',h.fs(end), ...
    'XGrid','on','YGrid','on');

set(ax(:,3),'TickDir','in');

for j=1:2
    ax(j,3).XTick = [1 200:200:ax(j,3).XTick(end)];
    ax(j,3).YTick = [1 200:200:ax(j,3).YTick(end)];
end

stfig_panlab(ax(:,1),{'A' 'B'},'hori','right','verti','bot','xoff',0,'yoff',0.075);
stfig_panlab(ax(:,1),{' standard DTW' ' DTW w/ symmetricP05 contraint'}, ...
    'hori','left','verti','bot','xoff',0,'yoff',0.075,'fontweight','normal');

stfig_panlab(ax(1,1:2),{'original signals','aligned signals'}, ...
    'hori','left','verti','top','xoff',0.01,'fontsize',h.fs(2),'fontweight','normal');

stfig_panlab(ax(1,3),{'distance matrix and warping curve'}, ...
    'hori','center','verti','bot','xoff',0,'fontsize',h.fs(2),'fontweight','normal','location','north');

xlabel(ax(end,2),'time index','fontsize',h.fs(3));
xlabel(ax(end,3),'{\itx} index','fontsize',h.fs(3));
ylabel(ax(end,3),'{\ity} index','fontsize',h.fs(3));

legend(ph(1,[1 2]),{'$x$','$y$'},'interpreter','latex','fontsize',h.fs(1),'location','southeast','NumColumns',2);
legend(pha(1,[1 2]),{'$\tilde{x}$','$\tilde{y}$'},'interpreter','latex','fontsize',h.fs(1),'location','southeast','NumColumns',2);

%%
h.printfig(mfilename);

end


