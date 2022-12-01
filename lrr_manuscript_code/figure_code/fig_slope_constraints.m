function [] = fig_slope_constraints()

dbstop if error; close all;
h = tw_helpers;

step_patterns = {'symmetric2','symmetricP05'};

S{1} = test_signals('slope_constraint_example');

D(1).X{1} = S{1}.X{1};
D(1).X{2} = S{1}.X{2};
D = repmat(D,2,1);

%
for i=1:length(D)
    [D(i).map,D(i).distance,D(i).map_info] = dtwm(D(i).X{1},D(i).X{2},'step_pattern',step_patterns{i});
    [D(i).Xa(1,:),D(i).Xa(2,:)] = align_signals(D(i).map,D(i).X{1},D(i).X{2},'aligntype','expansive');
end

%%
ax = stf([1 3 4; 2 3 4; 5 7 8; 6 7 8],[0.025 0.065 0.01 0.05],[0.075 0.075]);
ax = reshape(ax,[],2)';

for j=1:2
    shiftposy(ax(j,[1 2]),[-0.035 0.035]);
end

colors = lines(2);
lw = 3;
ls = {'-','-'};

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
    hd{i} = plot_distance(D(i).X{1},D(i).X{2},D(i).map);
    %set(hd{i}.mh,'LineWidth',3,'Color',pastelize([1 0 0],0.5));
    set(hd{i}.mh,'LineWidth',3);

    axes(ax(i,4));
    pat = select_step_pattern(step_patterns{i});
    plot_step_pattern(pat);

end

%%

axs = ax(:,[1 2]);
axis(axs,'tight');
axrescaley(0.05,axs);
set(axs,'XLim',[1 length(D(1).Xa)+1]);
set(axs,'XTick',[1 200:200:length(D(1).Xa)+1]);

set(ax(:,1),'XTickLabel',[]);
set(axs,'YTickLabel',[],'Box','off');

ax(1,4).XLim = ax(2,4).XLim;
ax(1,4).YLim = ax(2,4).YLim;

set(ax,'TickDir','out','ticklen',0.005*[1 1],'Fontsize',h.fs(end), ...
    'XGrid','on','YGrid','on');

set(ax(:,3),'TickDir','in');

for j=1:2
    ax(j,3).XTick = [1 ax(j,3).XTick];
    ax(j,3).YTick = [1 ax(j,3).YTick];
end

stfig_panlab(ax(:,1),{'A' 'B'},'hori','right','verti','bot','xoff',0,'yoff',0.075);
stfig_panlab(ax(:,1),{' standard DTW' ' DTW w/ local slope contraint'}, ...
    'hori','left','verti','bot','xoff',0,'yoff',0.075,'fontweight','normal');

stfig_panlab(ax(1,1:2),{'original signals','aligned signals'}, ...
    'hori','ce','verti','top','xoff',0.5,'fontsize',h.fs(3),'fontweight','normal');

stfig_panlab(ax(1,3),{'distance matrix and warping curve'}, ...
    'hori','center','verti','bot','xoff',0,'fontsize',h.fs(2),'fontweight','normal','location','north');

xlabel(ax(end,2),'time index','fontsize',h.fs(3));
xlabel(ax(end,3),'{\itx} index','fontsize',h.fs(3));
ylabel(ax(end,3),'{\ity} index','fontsize',h.fs(3));
xlabel(ax(end,4),'time index (relative)','fontsize',h.fs(3));

legend(ph(1,[1 2]),{'$x$','$y$'},'interpreter','latex','fontsize',h.fs(2),'location','southeast');
legend(pha(1,[1 2]),{'$\tilde{x}$','$\tilde{y}$'},'interpreter','latex','fontsize',h.fs(2),'location','southeast');

stfig_panlab(ax(1,4),ax(1,4).Title.String,'fontweight','normal', ...
    'verti','top','hori','left','xoff',0.01,'yoff',-0.05,'fontsize',h.fs(2));

stfig_panlab(ax(2,4),ax(2,4).Title.String,'fontweight','normal', ...
    'verti','top','hori','left','xoff',0.01,'yoff',-0.45,'fontsize',h.fs(2));

delete(ax(1,4).Title);
delete(ax(2,4).Title);

stfig_panlab(ax(1,4),{'step patterns'}, ...
    'hori','center','verti','bot','xoff',0,'fontsize',h.fs(2), ...
    'fontweight','normal','location','north');

%%
h.printfig(mfilename);

end


