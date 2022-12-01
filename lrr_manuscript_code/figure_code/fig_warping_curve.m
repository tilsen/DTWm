function [] = fig_warping_curve()

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
ax = stf([1 1 1 2 2 2; nan nan 3 3 3 3],[0.05 0.1 0.01 0.05],[0.10 0.15]);

max_len = 10;
tprops = {'Interpreter','latex','FontSize',h.fs(1),'hori','right'};

yo = 0.1;
xo = 0.5;
ms = 10;

for i=1:length(D)

    axes(ax(i));
    mm = D(i).map(:,1:min(max_len,end));

    H{i} = plot_time_map(mm,'label_mappings',true,'yoffset_text',0.2);

    set([H{i}.lh1 H{i}.lh2],'markersize',ms);
    
    set([H{i}.th],'Fontsize',h.fs(3));
    set([H{i}.th1 H{i}.th2],'Fontsize',h.fs(2));

    text(xo,0.5,'$k:$',tprops{:});
    text(xo,0-yo,'$j:$',tprops{:});
    text(xo,1+yo,'$i:$',tprops{:});

end

set([H{end}.th],'Fontsize',h.fs(2));

axis(ax,'tight');
axrescalex(0.05,ax);
axrescaley(0.35,ax);

for i=1:size(mm,2)
    ypos = [1 0]-2*yo*[-1 1];

    th(1,i) = text(mm(1,i),ypos(1),['$\phi_x(' num2str(i) ')$'],'verti','bot','hori','center', ...
        'interpreter','latex','fontsize',h.fs(2)+2);  

    th(2,i) = text(mm(2,i),ypos(2),['$\phi_y(' num2str(i) ')$'],'verti','top','hori','center', ...
        'interpreter','latex','fontsize',h.fs(2)+2);
end

th(1,2).String = {th(1,2).String,th(1,3).String,th(1,4).String};
delete(th(1,[3 4]));

th(2,7).String = {th(2,7).String,th(2,8).String,th(2,9).String};
delete(th(2,[8 9]));


text(-1,0.5,'$\phi(k)=\left(\phi_x(k),\phi_y(k)\right)$',tprops{:});

drawbrace([-0.95 1],[-0.95 0],0.01,'color','k');
set(gca,'Clipping','off')

%%

plh = stfig_panlab(ax,{'A','B','C'},'fontsize',h.fs(1),'xoff',[-0.05 -0.05 -0.1], ...
    'verti','mid','yoffset',[0 0 0.15]);


%%
h.printfig(mfilename);

end


