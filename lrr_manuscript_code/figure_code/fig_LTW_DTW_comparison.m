function [] = fig_LTW_DTW_comparison()

dbstop if error; close all
h = tw_helpers;

S{1} = test_signals('growth_decay_constant');
S{2} = test_signals('growth_decay_variable');

refix = 2;
for i=1:length(S)
    X = S{i}.X;
    for j=1:3
        A(i).map{j} = dtwm(X{refix},X{j});
        [~,A(i).dtw_aligned{j}] = align_signals(A(i).map{j},X{refix},X{j},'aligntype','reference');
    end
    [~,A(i).ltw_aligned] = ltwm(X,'length',length(X{refix}));
end

%%
ax = stf([4 2],[0.065 0.065 0.01 0.05],[0.025 0.05], ...
    'handlearray','matrix','aspect',1.1);

ls = {'-',':','-'};
lw = 3;
colors = lines(3);
colors = colors([1 2 2],:);
colors(2,:) = [0 0 0];

pargs = @(c){'linew',lw,'linesty',ls{c},'color',colors(c,:)};

oy = 0.03;

for i=1:2
    X = S{i}.X;
    rates = S{i}.rates;
    t = S{i}.t;
    dt = S{i}.dt;

    axes(ax(1,i));
    for j=1:length(X)
        pa = pargs(j);
        plot(t{j},X{j},pa{:}); hold on;
    end

    axes(ax(2,i));
    for j=1:length(X)
        pa = pargs(j);
        plot(t{j},rates{j}/S{i}.dt,pa{:}); hold on;
    end

    axes(ax(3,i));
    for j=1:3
        pa = pargs(j);
        plot(t{2},(j-2)*oy+A(i).ltw_aligned{j},pa{:}); hold on;
    end

    axes(ax(4,i));
    for j=1:3
        pa = pargs(j);
        plot(t{2},(j-2)*oy+A(i).dtw_aligned{j},pa{:}); hold on;
    end    

end

%%
axis(ax,'tight');

for i=1:2
    xlims = getlims(ax(:,i),'x');
    xlim(ax(:,i),xlims);
end

for i=1:size(ax,1)
    ylims = getlims(ax(i,:),'y');
    set(ax(i,:),'ylim',ylims);
end

axrescaley(0.075,ax);
axrescalex(0.025,ax);

set(ax,'Box','off','TickDir','out','ticklen',0.003*[1 1], ...
    'XGrid','on','YGrid','on','fontsize',h.fs(end));

set(ax(1:end-1,:),'XTickLabel',[]);
set(ax(:,2),'YTickLabel',[]);

ylabel(ax(1,1),'signal value (y)');
ylabel(ax(2,1),'rate (\Deltay/s)');
xlabel(ax(end,:),'time (s)');

obj_fontsize(gcf,'label',h.fs(3));

stfig_panlab(ax(1,:),{'constant-rate growth and decay' 'variable-rate growth and decay'}, ...
    'fontsize',h.fs(2)+4,'xoff',0,'hori','left','style','letter_title');

plabs = {'known process rates' 'LTW alignments' 'DTW alignments'};
stfig_panlab(ax(2:end,1),plabs,'xoff',0,'hori','left','fontsize',h.fs(3),'fontweight','normal');

%%
h.printfig(mfilename);

end


