function [] = tw_gen_gamdata()

dbstop if error;
h = tw_helpers;

t_margins = [0.1 0.1];

lrr_sets = {
    'CVCarl_deltas_100_P05'
    'CVCacl_deltas_100_P05'    };

for i=1:length(lrr_sets)
    load([h.lrr_dir 'lrr_' lrr_sets{i} '.mat']);

    X = []; D = [];

    G.su = grp2idx(G.subj);
    G.co = grp2idx(G.cond);

    for j=1:height(G)
        
        x = G.lrrs_gsd{j};
        d = [G.su(j) G.co(j)];

        for k=1:length(x)
            t = (0:(length(x{k})-1))/1000;
            X = [X; x{k}' t'];
            D = [D; k*ones(length(t),1) repmat(d,length(t),1)];
        end

    end

    D = uint8(D);
    X = single(X);

    R = table(X(:,1),X(:,2),D(:,1),D(:,2),D(:,3),'VariableNames',{'lrrgsd' 'time' 'trial' 'subj' 'cond'});

    trng = [min(R.time) max(R.time)];
    R = R(R.time>=(trng(1)+t_margins(1)) & R.time<=(trng(2)-t_margins(2)),:);

    writetable(R,[h.gam_dir 'gamdata_' lrr_sets{i} '.csv']);
    save([h.gam_dir 'gamdata_' lrr_sets{i} '.mat'],'R');
end


end
