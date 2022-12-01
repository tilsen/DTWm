function [] = fig_global_window()

dbstop if error; close all;
h = tw_helpers;

windows = {'','Itakura'};

S{1} = test_signals('global_window_example');

D(1).X{1} = S{1}.X{1};
D(1).X{2} = S{1}.X{2};
D = repmat(D,2,1);

%
for i=1:length(D)
    tic
    [D(i).map,D(i).distance,D(i).map_info] = dtwm(D(i).X{1},D(i).X{2},'window_type',windows{i});
    toc
    
    [D(i).Xa(1,:),D(i).Xa(2,:)] = align_signals(D(i).map,D(i).X{1},D(i).X{2},'aligntype','expansive');
end

%%
axpan = [1 1 1 3 3; 2 2 2 3 3; 4 4 4 6 6; 5 5 5 6 6];
ax = stf(axpan,[0.025 0.065 0.025 0.05],[0.075 0.075],'aspect',1.5);
ax = reshape(ax,[],2)';

for j=1:2
    shiftposy(ax(j,[1 2]),[-0.035 0.035]);
end

colors = lines(2);
lw = 2;
ls = {'-','--'};

yo = 0.1;

for i=1:length(D)

    axes(ax(i,1)); %#ok<*LAXES>
    for j=1:2
        ph(i,j) = plot(D(i).X{j}+yo*j,ls{j},'color',colors(j,:),'linew',lw); hold on;
    end

    axes(ax(i,2));
    for j=1:2
        pha(i,j) = plot(D(i).Xa(j,:)+yo*j,ls{j},'color',colors(j,:),'linew',lw); hold on;
    end

    axes(ax(i,3));
    %plot_distance(D(i).X{1},D(i).X{2},D(i).map);
    cost = D(i).map_info.localCostMatrix;
    cost(isnan(D(i).map_info.costMatrix)) = nan;
    stf_matrix(cost', ...
        'textvalues',false,'gridlines',false,'ydir','normal'); hold on;
    plot(D(i).map(1,:),D(i).map(2,:),'w-','linew',3);

end

%%

axs = ax(:,[1 2]);
axis(axs,'tight');
axrescaley(0.05,axs);
set(axs,'XLim',[1 length(D(1).Xa)+1]);
set(axs,'XTick',[1 200:200:length(D(1).Xa)+1]);

set(ax(:,1),'XTickLabel',[]);
set(axs,'YTickLabel',[],'Box','off');

set(ax,'TickDir','out','ticklen',0.003*[1 1],'Fontsize',h.fs(end), ...
    'XGrid','on','YGrid','on');

set(ax(:,3),'TickDir','in');

for j=1:2
    ax(j,3).XTick = [1 ax(j,3).XTick];
    ax(j,3).YTick = [1 ax(j,3).YTick];
end

stfig_panlab(ax(:,1),{'A' 'B'},'hori','right','verti','bot','xoff',0,'yoff',0.275);
stfig_panlab(ax(:,1),{'  standard DTW' '  DTW w/ global window constraint (Itakura)'}, ...
    'hori','left','verti','bot','xoff',0,'yoff',0.275,'fontweight','normal');

stfig_panlab(ax(1,1),{'original and aligned signals'}, ...
    'hori','left','verti','bot','xoff',0,'fontsize',h.fs(2),'fontweight','normal');

stfig_panlab(ax(1,3),{'distance matrix and warping curve'}, ...
    'hori','center','verti','bot','xoff',0,'fontsize',h.fs(2),'fontweight','normal','location','north');

xlabel(ax(end,2),'signal index','fontsize',h.fs(3));
xlabel(ax(end,3),'{\itx} index','fontsize',h.fs(3));
ylabel(ax(end,3),'{\ity} index','fontsize',h.fs(3));

legend(ph(2,[1 2]),{'$x$','$y$'},'interpreter','latex','fontsize',h.fs(2),'location','southeast');
legend(pha(2,[1 2]),{'$\tilde{x}$','$\tilde{y}$'},'interpreter','latex','fontsize',h.fs(2),'location','southeast');

%%
set(gcf,'InvertHardcopy','off','Color','w');
h.printfig(mfilename);

end


