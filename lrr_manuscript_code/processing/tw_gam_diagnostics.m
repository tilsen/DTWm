function [] = tw_gam_diagnostics(model)

dbstop if error; close all;

h = tw_helpers;
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

do_plot = true;

if nargin==0
    models = rdir([h.gam_dir 'gampred_*.mat']);
    models = regexp({models.name},'gampred_(\w+)\.mat$','tokens','once');
    for i=1:length(models)
        tw_gam_diagnostics(models{i}{:});
    end
    return;
else
    if contains(model,{'log' 'orig'})
        modeldata = model(1:max(regexp(model,'_'))-1);
    end
end

%%

kcheckfile = rdir([h.gam_dir 'gamkcheck_' model '.csv']);
x = readtable(kcheckfile.name);
x.name = repmat(model,height(x),1);
disp(x);

if ~do_plot, return; end

%%
file = rdir([h.gam_dir 'gampred_' model '.mat']);
load(file.name);


datafile = rdir([h.gam_dir 'gamdata_' modeldata '.mat']);
load(datafile.name);

Ds = Ds(~cellfun('isempty',Ds.subj),:);
Ds.cond = grp2idx(Ds.cond);
Ds.subj = grp2idx(Ds.subj);

R.trial = double(R.trial);
R.new_trial = [1; diff(R.trial)~=0];
ixs = find(R.new_trial);

R.pred = nan(height(R),1);
for i=1:length(ixs)
    pred = Ds.mu(Ds.subj==R.subj(ixs(i)) & Ds.cond==R.cond(ixs(i)),:);
    R.pred(ixs(i)+(0:length(pred)-1)) = pred;
end

if contains(model,'log')
    R.lrrgsd = log(R.lrrgsd);
    R.pred = log(R.pred);
end

R.resid = R.lrrgsd-R.pred;

R.resid_std = R.resid/std(R.resid);

%pd = makedist('tLocationScale','mu',0,'sigma',1,'nu',height(R)-1);

Nbins = 100;

[dens1,x1,y1] = histcounts2(R.pred,R.resid_std,Nbins);
[dens2,x2,y2] = histcounts2(R.pred,R.lrrgsd,Nbins);

%%
ax = stf([2 2],[0.05 0.075 0.01 0.05],[0.075 0.125],'aspect',1.15);

%qq-plot
axes(ax(1));
hh = qqplot(R.resid_std);
set(hh(1),'MarkerEdgeColor',0.5*ones(1,3),'Marker','.');
ylabel('residual quantiles');
axis square;
set(gca,'XTick',-5:5,'YTick',-5:5)

%density of residuals
axes(ax(3));
rr = minmax(R.resid_std');
rr = max(abs(rr))*[-1 1];
pnts = linspace(rr(1),rr(end),1000);

[f,xi] = ksdensity(R.resid_std,pnts);
plot(xi,f,'k-','linew',2);
title('density of residuals');
xlabel('standardized residual');
ylabel('density');
set(gca,'XTick',-5:5);

%residuals vs. linear pred
axes(ax(2));
imagesc(x1,y1,dens1'); hold on;
xlabel('prediction'); ylabel('standardized residual');
colormap('viridis');
plot(xlim,[0 0],'w--','linew',2);
title('residuals vs. prediction');

%response vs. fitted values
axes(ax(4));
imagesc(x2,y2,dens2'); hold on;
xlabel('prediction'); ylabel('empirical value');
colormap('viridis');
plot(xlim,xlim,'w--','linew',2);
title('empirical values vs. predictions');

%%
axis(ax,'tight');

set(ax(1),'ylim',getlimss(ax(1),'y'));

set(ax,'XGrid','on','YGrid','on','Fontsize',h.fs(end),'Box','off','ticklen',0.003*[1 1],'tickdir','out');

set(ax([2 4]),'YDir','normal');
set(ax([2 4]),'GridColor',0.95*ones(3,1));

set(ax(3),'YTickLabel',[]);

ax(1).Title.String = 'QQ plot of residuals vs. normal';

stfig_panlab(ax,{'i','iii','ii','iv'});

%%

h.printfig([mfilename '_' model]);

end




