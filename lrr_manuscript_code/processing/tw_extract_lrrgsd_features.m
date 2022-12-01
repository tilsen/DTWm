function [F] = tw_extract_lrrgsd_features(model)

dbstop if error; %close all;
h = tw_helpers;

load([h.gam_dir 'gampred_' model '.mat'],'Ds','Dd');

Ds.t = Ds.t-mean(Ds.t(1,:));
Dd.t = Dd.t-mean(Dd.t(1,:));

Fs = 1000;
sm_win = round(Fs*0.100);

%%
D = Dd;

D.onset1 = cellfun(@(c)c(1),D.cond1);
D.onset2 = cellfun(@(c)c(1),D.cond2);
D.coda1 = cellfun(@(c)c(end),D.cond1);
D.coda2 = cellfun(@(c)c(end),D.cond2);

%D = D(ismember(D.subj,'ALL'),:);

ff = [strcat('alpha',{'_max','_vel','_ons'}'); strcat('beta',{'_max','_vel','_ons'}')];

for i=1:length(ff)
    D.(ff{i}) = nan(height(D),1);
end

%% alpha feature
ixs = find(D.onset1=='0' & (D.onset2=='p'|D.onset2=='t') & D.coda1==D.coda2)';
D.ax(ixs) = true;

crit_alpha = stlm_selection_criteria('maxval',0.5,'danch',0.5);
crit_velmax = stlm_selection_criteria('maxval',1);

for i=ixs   
    x = smooth(D.mu(i,:),sm_win)';
    t = D.t(i,:);
    LM = stlm(x,t,'alpha_max',-0.150,'crit',crit_alpha,'ranget',[-inf 0]);
    LM(end+1) = stlm(x,t,'alpha_vel',LM(end).t,'crit',crit_velmax,'ranget',LM(end).t+[-0.200 0], ...
        'direction','backward','landmark','velmax');
    %LM(end+1) = stlm(x,t,'alpha_ons',LM(end).t,'direction','backward','landmark','velthr');

    D.alpha_max(i) = LM(1).t;    
    D.alpha_vel(i) = LM(2).t;    
    %D.alpha_ons(i) = LM(3).t;      
end

%% beta feature
ixs = find(D.coda1=='0' & (D.coda2=='p'|D.coda2=='t') & D.onset1==D.onset2)';
D.bx(ixs) = true;

crit_beta = stlm_selection_criteria('maxval',0.75,'danch',0.25);
crit_velmax = stlm_selection_criteria('maxval',0.5,'danch',0.5);

LM=[];
for i=ixs
    x = smooth(D.mu(i,:),sm_win)';
    t = D.t(i,:);

    LM = stlm(x,t,'beta',0.250,'crit',crit_beta,'ranget',[0 inf]);
    LM(end+1) = stlm(x,t,'beta_vel',LM(end).t,'crit',crit_velmax,'direction','backward','landmark','velmax');
    %LM(end+1) = stlm(x,t,'alpha_ons',LM(end).t,'direction','backward','landmark','velthr');

    D.beta_max(i) = LM(1).t;    
    D.beta_vel(i) = LM(2).t;    
    %D.beta_ons(i) = LM(3).t;  
end
F = D(D.ax | D.bx,:);



end

