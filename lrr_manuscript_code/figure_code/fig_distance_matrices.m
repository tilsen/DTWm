function [] = fig_distance_matrices()

dbstop if error; close all
h = tw_helpers;

S{1} = test_signals('cos_freq_var','params',{'n',3});
S{2} = test_signals('growth_decay_variable');
S{3} = test_signals('alignment_example');

D(1).X{1} = S{1}.X{1};
D(1).X{2} = S{1}.X{2};

D(2).X{1} = S{2}.X{1};
D(2).X{2} = S{2}.X{3};

D(3).X{1} = S{3}.X{1};
D(3).X{2} = S{3}.X{2};

for i=1:length(D)
    [D(i).map,D(i).distance,D(i).map_info] = dtwm(D(i).X{1},D(i).X{2});
end

%%
ax = stf([1 3 5; repmat([2 4 6],2,1)],[0.05 0.085 0.05 0.025],[0.05 0.095],'aspect',1.85);
ax = reshape(ax,2,[]);

params_short = {'signals_linespec',{'o-','s-'},...
'path_linespec','o-','matrixgrid',true};

params_long = {'signals_linespec',{'-','-'},...
'path_linespec','-','matrixgrid',false};

D(1).params = params_long;
D(2).params = params_long;
D(3).params = params_short;

D(1).ls = {'linestyle','-'};
D(2).ls = {'linestyle','-'};
D(3).ls = {'linestyle','-','markerfacecolor','w','marker','o'};

colors = lines(2);

for i=1:length(D)

    axes(ax(1,i));
    for j=1:2
        ph(i,j) = plot(D(i).X{j},'linew',2,D(i).ls{:},'color',colors(j,:)); hold on;
    end
    set(gca,'box','off','fontsize',h.fs(end),'YGrid','on','xgrid','on');

    axes(ax(2,i));
    H{i} = plot_dtw_matrix(D(i).map_info,D(i).X{1},D(i).X{2},'parent',gca,D(i).params{:}); hold on;

    set([H{i}.s1 H{i}.s2 H{i}.maph],'linew',2);

    if i<length(D)
        delete(H{i}.cbh);
    else
        shiftposx(H{i}.cbh,-0.015);
        H{i}.cbh.Position(3) = 0.01;
        set(H{i}.cbh,'Fontsize',h.fs(end));
    end

    set(H{i}.axs,'Fontsize',h.fs(end));
    set(H{i}.axs(1),'YTickLabel',[]);
    switch(i)
        case {1,2}
            H{i}.axs(1).XTick = unique([1 H{i}.axs(1).XTick]);
            H{i}.axs(2).YTick = unique([1 H{i}.axs(2).YTick]);
            
        otherwise
            H{i}.axs(1).XTick = 1:max(H{i}.axs(1).XLim);
            set(H{i}.s1,'markerfacecolor','w','marker','o');
            set(H{i}.s2,'markerfacecolor','w','marker','o');
            set(H{i}.maph,'markerfacecolor','k','marker','o');
            text(1,2,'$d_{1,2}$','Interpreter','latex', ...
                'FontSize',h.fs(end),'color','w','hori','center', ...
                'parent',H{i}.axd,'fontweight','bold');
    end
    set(H{i}.axs(2),'XTickLabel',[]);
    if i>1
       delete(H{i}.ylabs2);
    end

end

%%

stfig_panlab(ax(1,:),{'A','B','C'},'fontsize',h.fs(1),'xoff',-0.1,'verti','mid','yoffset',0);

axis(ax(1,:),'tight');
xlabel(ax(1,:),'sample index','FontSize',h.fs(end));

for i=1:2
    ax(1,i).XTick = [1 ax(1,i).XTick];
end

legh = legend(ph(1,:),{'$x$','$y$'},'location','southeast','fontsize',h.fs(2),'Interpreter','latex');
shiftposx(legh,0.025);

%%
h.printfig(mfilename);

end


