function [] = fig_alignment_methods()

dbstop if error; close all
h = tw_helpers;

S{1} = test_signals('alignment_example');

D(1).X{1} = S{1}.X{1};
D(1).X{2} = S{1}.X{2};
[D(1).map,D(1).distance,D(1).map_info] = dtwm(D(1).X{1},D(1).X{2});

D = repmat(D,4,1);

map = D(1).map;

% compression
D(2).map = map(:,~any(diff([map [0;0]],[],2)==0));

%asymmetric
ix1 = unique(map(1,:));
ix2 = arrayfun(@(c)map(2,find(map(1,:)==c,1,'last')),ix1);
D(3).map = [ix1; ix2];

ix2 = unique(map(2,:));
ix1 = arrayfun(@(c)map(1,find(map(2,:)==c,1,'last')),ix2);
D(4).map = [ix1; ix2];

%
for i=1:length(D)
    D(i).Xa(1,:) = D(i).X{1}(D(i).map(1,:));
    D(i).Xa(2,:) = D(i).X{2}(D(i).map(2,:));
end

%%
axo = stf([2 2],[0.025 0.065 0.01 0.05],[0.075 0.10]);

for j=1:length(axo)
    ax(j,:) = stfig_subaxpos(axo(j),[1 1 2 2 2 3 3]',[0 0 0 0 0.01 0.01]);
end
delete(axo);

ax = ax([1 3 2 4],:);

ms = 6;
colors = lines(2);
ls = {'o-','s:'};

for i=1:length(D)

    axes(ax(i,1));
    for j=1:2
        ph(i,j) = plot(D(i).X{j},ls{j},'color',colors(j,:),'linew',3,'markerfacecolor','w'); hold on;
    end
    
    axes(ax(i,3));
    for j=1:2
        pha(i,j) = plot(D(i).Xa(j,:),ls{j},'color',colors(j,:),'linew',3); hold on;
    end
    
    axes(ax(i,2));
    Hx{i} = plot_time_map(D(1).map,'markers',{'o','s'},'label_indices',false); hold on;
    set([Hx{i}.lh1 Hx{i}.lh2],'markerfacecolor',[.95 .95 .95],'Color',[.5 .5 .5]);
    set(Hx{i}.ch,'linestyle',':');

    H{i} = plot_time_map(D(i).map,'label_mappings',true,'markers',{'o','s'});

    set([H{i}.lh1 H{i}.lh2],'markersize',ms);
    set([H{i}.th H{i}.th1 H{i}.th2],'Fontsize',h.fs(end));
    delete([H{i}.th1 H{i}.th2]);

end

%%
set(ax,'Box','off','TickDir','in','ticklen',0.005*[1 1],'Fontsize',h.fs(end));

axis(ax(:,2),'tight');
axrescaley([-0.2 0.4],ax(:,2));

xlim(ax,[1 10]);

set(ax(:,[1 3]),'XGrid','on','YTick',[],'YColor','none');

for i=1:size(ax,1)
    L = max([length(D(i).X{1}) length(D(i).X{2}) size(D(i).Xa,2)]);
    set(ax(i,:),'XTick',1:L);
end

set(ax(1,1),'XTick',1:8)

legend(ph(1,:),{'$x$','$y$'},'location','northeast','fontsize',h.fs(2),'interpreter','latex');
legend(pha(1,:),{'$\tilde{x}$','$\tilde{y}$'},'location','northwest','fontsize',h.fs(2),'interpreter','latex');

plh = stfig_panlab(ax(:,1),{'A' 'B' 'C' 'C^{\prime}'},'fontsize',h.fs(1),'xoff',-0.05, ...
    'verti','bot','yoffset',0.01,'hori','left','fontweight','normal');

strs = {'expansive (external observer)','compressive (external observer)','x''s reference frame','y''s reference frame'};
for i=1:length(plh)
    plh(i).String = ['{\bf' plh(i).String '}' blanks(2) strs{i}];
end

xlabel(ax(2,3),'time index','fontsize',h.fs(3));

%%
h.printfig(mfilename);

end


