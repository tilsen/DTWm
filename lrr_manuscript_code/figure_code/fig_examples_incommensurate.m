function [] = fig_examples_incommensurate()

dbstop if error; close all
h = tw_helpers;

S{1} = test_signals('cos_per_incom');
S{2} = test_signals('growth_incom');
S{3} = test_signals('similar_shape');

%%
ax = stf([1 3],[0.035 0.05 0.01 0.05],[0.05 0],'aspect',3);

for i=1:length(S)

    colors = 0.5*ones(length(S{i}.X),3);
    colors(end,:) = [1 0 0];
    lw = 2*ones(length(S{i}.X),1);
    lw(end) = 3;

    axes(ax(i));
    X = S{i}.X;
    t = S{i}.t;
    for j=1:length(X)
        plot(t{j},X{j},'linew',lw(j),'color',colors(j,:)); hold on;
    end
end

%%
axis(ax,'tight');
axrescalex(0.025,ax);
axrescaley(0.05,ax);

set(ax,'XGrid','on','YGrid','on','Tickdir','out',...
    'ticklen',0.003*[1 1],'fontsize',h.fs(end),'box','off');

set(ax,'XTicklabel',[],'YTickLabel',[]);

stfig_panlab(ax,{'A','B','C'},'hori','right','yoff',0,'xoff',-0.025,'verti','mid','fontsize',h.fs(1));

%%
h.printfig(mfilename);

end

