function [] = fig_LTW()

dbstop if error; close all;
h = tw_helpers;

L_orig = 11;
dt = 1/(L_orig-1);

t = 0:dt:1;
x = sin(2*pi*t);

L_warp = [15 7];

Wfac = L_warp/L_orig;

for i=1:length(L_warp)
    ltw(i).ti = linspace(t(1),t(end),L_warp(i));
    ltw(i).tn = t(1):dt:(L_warp(i)-1)*dt;
    ltw(i).xi = interp1(t,x,ltw(i).ti);
    ltw(i).w = (L_warp(i)-1)/(L_orig-1);
end

%%
ax = stf([1 1 2 3 3 4 4; 5 5 6 7 7 8 8]',[0.075 0.075 0.01 0.05],[0.075 0.00]);
ax = reshape(ax,[],2);

colors = [0 0 0; lines(2)];

txo = -0.075;
tyo = -0.1;

ms = 8;

for i=1:2

    ti = ltw(i).ti;
    xi = ltw(i).xi;
    tn = ltw(i).tn;

    axes(ax(1,i));
    pho(i) = plot(t,x,'o-','color',colors(1,:), ...
        'markerfacecolor','w','linew',2,'markersize',ms); hold on;

    text(2*txo,0,'$x(t)$','Interpreter','latex','HorizontalAlignment','right','FontSize',h.fs(2));
    

    axes(ax(2,i));
    for j=1:length(t)
        plot(t(j),0,'o','Color','k', ...
            'markerfacecolor','w','markersize',ms); hold on;
        text(t(j),tyo,num2str(t(j)),'verti','top','hori','center', ...
            'fontsize',h.fs(end));
    end

    text(txo,0,'$t_{orig}$','Interpreter','latex','HorizontalAlignment','right','FontSize',h.fs(1));

    axes(ax(3,i));
    
    if i==1
        tiprops = {'hori','left','rot',45,'verti','bot'};
    else
        tiprops = {'hori','center','rot',0,'verti','bot'};
    end

    for j=1:length(ti)
        line([ti(j) tn(j)],[0 -1],'color','k','linestyl','--','linew',1); hold on;
    end

    for j=1:length(ti)
        plot(ti(j),0,'o','Color','k', ...
            'markerfacecolor','w','markersize',ms); hold on;
        text(ti(j),-tyo,[num2str(ti(j),'%1.2f')],tiprops{:}, 'fontsize',h.fs(end));   

        plot(tn(j),-1,'o','Color','k', ...
            'markerfacecolor',colors(i+1,:),'markersize',ms);
        text(tn(j),-1+tyo,num2str(tn(j),'%1.2f'),'verti','top','hori','center','fontsize',h.fs(end));            
    end

    text(txo,0,'$t_{interp}$','Interpreter','latex','HorizontalAlignment','right','FontSize',h.fs(1));
    text(txo,-1,'$\tilde{t}_{norm}$','Interpreter','latex','HorizontalAlignment','right','FontSize',h.fs(1));

    axes(ax(4,i));
    ph(1,i) = plot(t,x,'o-','color',colors(1,:), ...
        'markerfacecolor','w','markersize',ms); hold on;
    ph(2,i) = plot(tn,xi,'o-','color','k', ...
        'markerfacecolor',colors(i+1,:),'markersize',ms); hold on;

    %text(txo,0,'$\tilde{x}_{warped}$','Interpreter','latex','HorizontalAlignment','right','FontSize',h.fs(3));    

end

%%

axis(ax,'tight');
set(ax,'XLim',[-0.075 1.45],'Box','off','fontsize',h.fs(end),'tickdir','out','ticklen',0.003*[1 1]);

set(ax(3,:),'YLim',[-2 0.5]);

set(ax([2 3],:),'Visible','off');

set(ax(1,:),'XColor','none','XGrid','on','YTick',-1:0.5:1,'YGrid','on','YLim',[-1.1 1.1]);

set(ax(end,:),'XGrid','on','YGrid','on');
axrescaley(0.05,ax(end,:));

xlabel(ax(end,:),'time (s)','fontsize',h.fs(3));

for i=1:2
    lhx(i) = legend(ph(:,i),{'$x_{orig}$','$\tilde{x}_{warped}$'}, ...
        'location','northeast','fontsize',h.fs(1),'Interpreter','latex');

    lhx(i).Position(2) = lhx(i).Position(2)+0.025;

    legend(pho(i),{'$x_{orig}$'}, ...
        'location','northeast','fontsize',h.fs(1),'Interpreter','latex');    

    text(-0.05,-0.75,['$\mathcal{W}=' num2str(ltw(i).w,'%1.1f') '$'],'parent',ax(end,i),'Interpreter','latex','FontSize',h.fs(1));
end

stfig_panlab(ax(1,:),{'A  linear stretching','B  linear compression'},'fontsize',h.fs(1),'xoff',-0.1,'hori','left');

%%
h.printfig(mfilename);

end


