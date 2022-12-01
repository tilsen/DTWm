function [] = fig_temporal_structure()

dbstop if error; close all;
h = tw_helpers;

winsize = 51; %local slope window size
params = {'n',8};

mean_type = 'geometric'; sdf = 'gsd'; muf = 'gm';
%mean_type = 'harmonic'; sdf = 'hsd'; muf = 'hm';

S = test_signals('growth_decay_phases',params);
X = table(S.X',S.rates',S.t','VariableNames',{'X','rates','t'});

overwrite = false;

ixs{1} = 1:(height(X)/2);
ixs{2} = (height(X)/2+1):height(X);

geom_std = @(x)exp(std(log(x)));

%%
datafile = ['.' filesep 'figure_code' filesep 'lrr_example_complex.mat'];
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
axpan = [repmat([1 2],6,1); repmat([3 4],4,1); repmat([5 6],4,1)];
axpan(6,1) = 3;

ax = stf(axpan,[0.065 0.065 0.01 0.035],[0.15 0.10],'aspect',0.85);

colors = hsv(25);
colors = colors([1:4 floor(end/2)+(1:4)],:);
lw = 2;
ls = reshape(repmat({'-' '-'},4,1),[],1);

%% orignal signals
axes(ax(1));
for j=1:Ns
    ph(j) = plot(X{j},ls{j},'color',colors(j,:),'linew',lw); hold on;
end

for j=1:Ns
    [~,locs(j)] = findpeaks(X{j});
end
locs = fliplr(reshape(locs([1 4 5 8]),2,[])');
sets = {'A' 'B'};
for i=1:2
    drawbrace([locs(i,1) 1.01],[locs(i,2) 1.01],0.0035,'color','k','linew',2);
    text(mean(locs(i,:)),1.01,sets{i},'fontsize',h.fs(3),'hori','center','verti','bot');
end


%% matrix of LRRs
axes(ax(2));
[phmm, axm] = plot_lrr_array(LRR,colors);

%% process rates
axes(ax(3));
for i=1:Ns
    plot(S.rates{i}/S.dt,ls{i},'color',colors(i,:),'linew',lw); hold on;
end

%% mean of linearly normalized LRRs
axes(ax(4));
tix = linspace(0,1,length(N.gsd));
phm(1) = plot(tix,N.(muf),'-','color','k','linew',lw); hold on;

cix = [1 size(colors,1)];
for i=1:2
    cond_LRRn = T.LRRn(:,ixs{i},:);
    cond_LRRn = reshape(cond_LRRn,[],size(cond_LRRn,3));
    cond_mu{i} =  geomean(cond_LRRn);
    cond_sd{i} =  geom_std(cond_LRRn);
    phm(i+1) = plot(tix,cond_mu{i},'color',colors(cix(i),:),'linew',lw);
end    

%% deviation
axes(ax(5));
for i=1:Ns
    plot(M.(sdf){i},ls{i},'color',colors(i,:),'linew',lw); hold on;
end

%% sd of linearly normalized LRRs
axes(ax(6));
phs(1) = plot(tix,N.(sdf),'-','color','k','linew',lw); hold on;
for i=1:2
    phs(i+1) = plot(tix,cond_sd{i},'color',colors(cix(i),:),'linew',lw);
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

%set(ax([3 4]),'ylim',getlims(ax([3 4]),'y'));
set(ax([5 6]),'ylim',getlims(ax([5 6]),'y'));

axrescale(ax(tsix),0.025,0.05);
axrescale(ax(1),[],0.01);

set(ax(2),'Visible','off');
set([ax(2) axm(1,:)],'XAxisLocation','top');
set([ax(2) axm(1,:)],'XAxisLocation','top');
for i=1:size(axm,1)
    ylabel(axm(i,1),sprintf('%1.0f ',i), ...
        'rotation',0,'hori','right','fontsize',h.fs(end),'verti','mid');
    xlabel(axm(1,i),sprintf('%1.0f',i), ...
        'rotation',0,'hori','center','fontsize',h.fs(end),'verti','bot');    
end
xlabh = xlabel(ax(2),'reference signal','fontsize',h.fs(3),'Visible','on');
ylabh = ylabel(ax(2),'query signal','fontsize',h.fs(3),'Visible','on');

xlabh.Position(2) = xlabh.Position(2)-0.035;
ylabh.Position(1) = ylabh.Position(1)+0.035;

xlabel(ax(5),'signal time index','FontSize',h.fs(end));
xlabel(ax(6),'normalized event time','FontSize',h.fs(end));

ylabel(ax(1),'arbitrary state value','fontsize',h.fs(end));

legend(phm,{'all signals' 'subset A' 'subset B'}, ...
    'fontsize',h.fs(end)-2,'location','northeast');

yo = 0.04;
xo = -0.05;
stfig_panlab(ax,{'A' 'B' 'C' 'D' 'E' 'F'}, ...
    'fontsize',h.fs(2),'xoff',xo,'yoff',yo,'hori','right','verti','bot');

labs = {'signals','',...
    'generating process rates',...
    ['LRR ' mean_type ' mean'],...
    ['LRR ' mean_type ' stdev (by signal)'],...
    ['LRR ' mean_type ' stdev']};

stfig_panlab(ax,labs, ...
    'fontsize',h.fs(3),'fontweight','normal','hori','left','xoff',0.025,'yoff',yo,'verti','bot');

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
