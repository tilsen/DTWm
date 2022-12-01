function [] = tw_lrr_batch()

dbstop if error;
h = tw_helpers;

h.overwrite = true;

h.use_parallel = true;
h.poolsize = 12; %number of cores to use

h.windur = 0.100;
h.dt = 1e-3;
h.winsize = round(h.windur/h.dt);


%% analyses
A(1).lrr_set = 'CVCarl_deltas_100_P05';
A(end).dtw_params = {'step_pattern','symmetricP05'};

A(end+1).lrr_set = 'CVCacl_deltas_100_P05';
A(end).dtw_params = {'step_pattern','symmetricP05'};

A = struct2table(A);
A.use_deltas = true(height(A),1);
A.reverse_time = false(height(A),1);

%%
for i=1:height(A)
    process_localrate(A(i,:),h);
end

end

%%
function [] = process_localrate(A,h)

use_parallel = h.use_parallel;
poolsize = h.poolsize;

lrr_set = A.lrr_set{1};

pp = strsplit(lrr_set,'_');
dataset = pp{1};

if ~exist([h.lrr_dir lrr_set],'dir'), mkdir([h.lrr_dir lrr_set]); end

load([h.datasets_dir 'dataset_' dataset '.mat'],'D');

G = D.Properties.UserData{1};
IXfcn = D.Properties.UserData{2};

PAR.lrr_set = lrr_set;
PAR.winsize = h.winsize;
PAR.dtw_params = A.dtw_params(1,:);
PAR.use_deltas = A.use_deltas;
PAR.reverse_time = A.reverse_time;

ixs = arrayfun(@(c){IXfcn(c)},G.id);

if ~iscell(D.X)
    D.X = arrayfun(@(c){squeeze(D.X(c,:,:))},(1:height(D))');
end

G.X = cellfun(@(c){D.X(c)},ixs);

if PAR.use_deltas
    G.X = cellfun(@(c){cellfun(@(d){row_gradients(d)},c)},G.X);
end

%reverse times (for use with open-end)
if PAR.reverse_time
    G.X = cellfun(@(c){cellfun(@(d){fliplr(d)},c)},G.X); %#ok<*UNRCH>
end

%convert to double
try
    if isa(G.X{1}{1},'single')
        G.X = cellfun(@(c){cellfun(@(d){double(d)},c)},G.X);
    end
catch
end

sstr = 'processing local relative rates: %s\t%03i/%03i';

switch(use_parallel)    
    case true
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj), parpool(poolsize); end
end

tic
for i=1:height(G)
    status_str = status(sprintf(sstr,lrr_set,i,height(G))); %#ok<NASGU> 

    outfile = [h.lrr_dir lrr_set filesep sprintf('lrr_%03i.mat', G.id(i))];
    if ~h.overwrite && exist(outfile,'file')
        continue;
    end

    %lrr target signals relative to comparisons:
    switch(use_parallel)
        case false
            [loc_rel_rates,distances] = lrr_pairwise(G.X{i},'winsize',PAR.winsize,PAR.dtw_params{:});
        case true
            [loc_rel_rates,distances] = lrr_pairwise(G.X{i},'useparallel',true,'winsize',PAR.winsize,PAR.dtw_params{:});
    end
    g = G(i,:); g.X = [];
    save(outfile,'loc_rel_rates','g','distances','PAR');
end
elapsed = toc;
fprintf('\t elapsed time: %1.1f seconds',elapsed);
status('reset');

end

%%
function [dX] = row_gradients(X)

dX = nan(size(X));
for i=1:size(X,1)
    dX(i,:) = gradient(X(i,:));
end

end
