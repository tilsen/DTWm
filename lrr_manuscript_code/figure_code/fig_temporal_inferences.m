function [] = fig_temporal_inferences()

dbstop if error; close all;
h = tw_helpers;

winsize = 101; %local slope window size
params = {'n',7};

mean_type = 'geometric'; sdf = 'gsd'; muf = 'gm';
%mean_type = 'harmonic'; sdf = 'hsd'; muf = 'hm';

S = test_signals('growth_set',params);
X = table(S.X',S.rates',S.t','VariableNames',{'X','rates','t'});

Ns = height(X);

overwrite = false;

%%
datafile = ['.' filesep 'figure_code' filesep 'lrr_example.mat'];
if overwrite || ~exist(datafile,'file')
    T = lrr_analysis(X,winsize);
    save(datafile,'T');
else
    load(datafile,'T');
end

T.X = X.X;
ff = fieldnames(T);
for i=1:length(ff)
    eval([ff{i} '=T.' ff{i} ';'])
end

%%
axpan = [repmat([1 2],3,1); repmat([3 4],2,1); repmat([5 6],2,1)];
ax = stf(axpan,[0.065 0.065 0.01 0.075],[0.15 0.10],'aspect',0.85);

colors = hsv(Ns);
lw = 2;
ls = '-';

%% orignal signals
axes(ax(1));
for j=1:Ns
    ph(j) = plot(X{j},ls,'color',colors(j,:),'linew',lw); hold on;
    ths(j) = text(length(X{j}),X{j}(end),[' ' num2str(S.constant_rates(j),'%1.2f')], ...
        'hori','left','verti','mid','fontsize',h.fs(end),'rotation',90);
end

%% matrix of LRRs
axes(ax(2));
[phmm, axm] = plot_lrr_array(LRR,colors);

%% means/sds by signal
for i=1:Ns
    plot(M.(muf){i},ls,'color',colors(i,:),'linew',lw,'parent',ax(3)); hold(ax(3),'on');
    plot(M.(sdf){i},ls,'color',colors(i,:),'linew',lw,'parent',ax(5)); hold(ax(5),'on');
end

%% sd/means of linearly normalized LRRs
muff = {'gm' 'am' 'hm'};
sdff = {'gsd' 'asd' 'hsd'};
mulabs = {'geometric' 'arithmetic' 'harmonic'};
colorsm = [0 0 0; .5 .5 .5;  1 0 0];
lsm = {'-' '-' '--'};

tix = linspace(0,1,length(N.(muff{1})));

for i=1:length(muff)
    phm(i) = plot(tix,N.(muff{i}),lsm{i},'color',colorsm(i,:),'linew',2,'parent',ax(4)); hold(ax(4),'on');
    phs(i) = plot(tix,N.(sdff{i}),lsm{i},'color',colorsm(i,:),'linew',2,'parent',ax(6)); hold(ax(6),'on');
end

%%
set(ax,'Box','off',h.ticks{:},'fontsize',h.fs(end)-2);

tsix = [1 3 4 5 6];
set(ax(1),'YTickLabel',[]);
set(ax(tsix),'XGrid','on','YGrid','on');
axis(ax(tsix),'tight'); 

axts = ax([1 3 5]);
axtsn = ax([4 6]);

set(axts,'xlim',getlims(axts,'x'));
set(axtsn,'xlim',getlims(axtsn,'x'));

set(ax([3 4]),'ylim',getlims(ax([3 4]),'y'));
set(ax([5 6]),'ylim',getlims(ax([5 6]),'y'));

axrescale(ax(tsix),0.025,0.05);
axrescale(ax(1),[],[-0.01 0.10]);

set([ax(2) axm(1,:)],'XAxisLocation','top');
for i=1:size(axm,1)
    ylabel(axm(i,1),sprintf('%1.2f ',S.constant_rates(i)), ...
        'rotation',0,'hori','right','fontsize',h.fs(end),'verti','mid');
    xlabel(axm(1,i),sprintf('%1.2f',S.constant_rates(i)), ...
        'rotation',0,'hori','center','fontsize',h.fs(end),'verti','bot');    
end

set(ax(2),'Visible','off');
xlabh = xlabel(ax(2),'reference signal rate','fontsize',h.fs(3),'Visible','on');
ylabh = ylabel(ax(2),'query signal rate','fontsize',h.fs(3),'Visible','on');

xlabh.Position(2) = xlabh.Position(2)-0.025;
ylabh.Position(1) = ylabh.Position(1)-0.025;

legh = legend(phm,mulabs, ...
    'fontsize',h.fs(end),'location','northeast','NumColumns',1,'AutoUpdate','off');

legh.Position(2)=legh.Position(2)+0.025;

xlabel(ax([5]),'signal time index','FontSize',h.fs(end));
xlabel(ax([6]),'normalized event time','FontSize',h.fs(end));

%ylabel(ax(1),'arbitrary state value','fontsize',h.fs(end));

yo = 0.04;
xo = -0.05;
stfig_panlab(ax,{'A' 'B' 'C' 'D' 'E' 'F'}, ...
    'fontsize',h.fs(2),'xoff',xo,'yoff',yo,'hori','right');

labs = {'signals','',...
    ['LRR ' mean_type ' mean (by signal)'],...
    'LRR mean',...
    ['LRR ' mean_type ' stdev (by signal)'],...
    ['LRR stdev']};
stfig_panlab(ax,labs, ...
    'fontsize',h.fs(3),'fontweight','normal','hori','left','xoff',0.025,'yoff',yo);

%%
h.printfig(mfilename);

end

%%
function [T] = lrr_analysis(S,winsize)

X = S.X;
rates = S.rates;
Ns = length(X);

for a=1:Ns
    for b=1:Ns
        [D(a,b).map,D(a,b).dist,~] = dtwm(X{a},X{b});
        [D(a,b).Xa(1,:),D(a,b).Xa(2,:),D(a,b).map_ref] = align_signals(D(a,b).map,X{a},X{b},'aligntype','reference');
        D(a,b).LRR = lrr(D(a,b).map,'winsize',winsize);
        D(a,b).relrate = rates{a}./rates{b}(D(a,b).map_ref(2,:));
    end
end


LRR = reshape({D.LRR},Ns,[])';

L = cellfun('length',LRR);
maxL = max(L(:));
medL = median(L(:));

for a=1:size(LRR,1)
    for b=1:size(LRR,2)
        LRRm(a,b,:) = [LRR{a,b} nan(1,maxL-L(a,b))];
    end
end

for a=1:size(LRR,1)
    for b=1:size(LRR,2)
        [~,x] = ltwm(LRR{a,b},'length',medL);
        LRRn(a,b,:) = [x{:}];
    end
end

T.Ns = Ns;
T.X = X;
T.LRR = LRR;
T.LRRm = LRRm;
T.LRRn = LRRn;

T.L = L;
T.medL = medL;
T.maxL = maxL;
T.minLRR = min([LRR{:}]);
T.maxLRR = max([LRR{:}]);
T.rangeLRR = T.maxLRR;

harm_mean = @(x)size(x,1)./sum(1./x);
harm_std = @(x)sqrt(    sum(((1./x)-(1./harm_mean(x))).^2)  /  (size(x,1))    );

geom_mean = @(x)prod(x).^(1/size(x,1));
%geom_std = @(x)exp(std(log(x)));
geom_std = @(x)exp( sqrt(   sum(    (log(x./geom_mean(x)).^2) / size(x,1)) )     );

M.gm = arrayfun(@(c){geom_mean(cell2mat(LRR(:,c)))},(1:Ns));
M.gsd = arrayfun(@(c){geom_std(cell2mat(LRR(:,c)))},(1:Ns));
N.gm = geom_mean(reshape(LRRn,Ns*Ns,[]));
N.gsd = geom_std(reshape(LRRn,Ns*Ns,[]));

M.hm = arrayfun(@(c){harm_mean(cell2mat(LRR(:,c)))},(1:Ns));
M.hsd = arrayfun(@(c){harm_std(cell2mat(LRR(:,c)))},(1:Ns));
N.hm = harm_mean(reshape(LRRn,Ns*Ns,[]));
N.hsd = harm_std(reshape(LRRn,Ns*Ns,[]));

M.am = arrayfun(@(c){mean(cell2mat(LRR(:,c)))},(1:Ns));
M.asd = arrayfun(@(c){std(cell2mat(LRR(:,c)))},(1:Ns));
N.am = mean(reshape(LRRn,Ns*Ns,[]));
N.asd = std(reshape(LRRn,Ns*Ns,[]));

T.M = M;
T.N = N;

end
