function [] = tw_summarize_traj()

dbstop if error;
h = tw_helpers;

dataset = 'CVCarl';

load([h.datasets_dir 'dataset_' dataset '.mat']);

subjs = unique(D.subj);
conds = unique(D.cond);

Dk = [];
for i=1:length(subjs)
    for j=1:length(conds)
        Dx = tabindex(D,'subj',subjs{i},'cond',conds{j});

        Dk(end+1).subj = subjs{i};
        Dk(end).cond = conds{j};
        Dk(end).mu = nanmean(Dx.X,1);
        Dk(end).sd = nanstd(Dx.X,[],1);
        Dk(end).se = Dk(end).sd./sum(~isnan(Dx.X));
        Dk(end).p95 = prctile(Dx.X,[2.5 97.5]);
        Dk(end).t = Dx.Properties.UserData{3}.t_kin;

    end
end

%all subj
for j=1:length(conds)
    Dx = tabindex(D,'cond',conds{j});

    Dk(end+1).subj = 'ALL';
    Dk(end).cond = conds{j};
    Dk(end).mu = nanmean(Dx.X,1); %#ok<*NANMEAN> 
    Dk(end).sd = nanstd(Dx.X,[],1); %#ok<*NANSTD> 
    Dk(end).se = Dk(end).sd./sum(~isnan(Dx.X));
    Dk(end).p95 = prctile(Dx.X,[2.5 97.5]);
    Dk(end).t = Dx.Properties.UserData{3}.t_kin;

end

D = struct2table(Dk);
D.Properties.UserData{1} = Dx.Properties.UserData{3};

%%
save([h.datasets_dir 'summary_' dataset '.mat'],'D');

end