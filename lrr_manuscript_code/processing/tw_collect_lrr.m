function [] = tw_collect_lrr(varargin)

dbstop if error; close all;
h = tw_helpers;

%threshold for outlier exclusion for statistical summaries:
zthresh = norminv(0.99);

overwrite = true;

if nargin==0
    lrr_sets = {
        'CVCarl_deltas_100_P05' 
        'CVCacl_deltas_100_P05'};      

    for i=1:length(lrr_sets)
        tw_collect_lrr(lrr_sets{i});
    end
    return;

elseif nargin==1
    lrr_set = varargin{1};
end


outfile = [h.lrr_dir 'lrr_' lrr_set '.mat'];
if ~overwrite && exist(outfile,'file'), return; end

%not numerically stable:
%geom_mean = @(x)prod(x).^(1/size(x,1));

%numerically stable:
%geom_mean = @(x)exp( sum(log(x)) / size(x,1));
%geom_std = @(x)exp( sqrt(   sum(    (log(x./geom_mean(x)).^2) / size(x,1)) )     );

%numerically stable and handles nan:
geom_mean = @(x)exp( nansum(log(x)) ./ sum(~isnan(x)) ); %#ok<*NANSUM> 
geom_std = @(x)exp( sqrt(   nansum((log(x)-log(geom_mean(x))).^2) ./ sum(~isnan(x))  )     );


%%

pp = strsplit(lrr_set,'_');
dataset = pp{1};

datasetfile = [h.datasets_dir 'dataset_' dataset '.mat'];
if ~exist(datasetfile,'file'), fprintf('dataset: %s not found\n',datasetfile); return; end
load(datasetfile,'D');

%local rate calculation groups
G = D.Properties.UserData{1};
IXfcn = D.Properties.UserData{2};
Gi = D.Properties.UserData{3};

%%
files = rdir([h.lrr_dir lrr_set filesep 'lrr_*.mat']);

if isempty(files), fprintf('dataset %s: lrr files not found\n',datasetfile); return; end

for i=1:length(files)


    status_str = status(sprintf('%s: %i/%i',lrr_set,i,length(files)));
    x = load(files(i).name);

    %array of pairwise loc_slopes: slopes of column signal (query/comparison) to row signal (ref/target)
    loc_rel_rates = x.loc_rel_rates;

    %transpose ref (target) / query (comparison) so that rows are target signals
    LRR = loc_rel_rates';

    %correct too low/negative slopes (only necessary in the absence of a
    %local slope constraint):
    maxrate = max([LRR{:}]);
    minrate = 1/(2*maxrate);

    LRRcorrected = zeros(size(LRR));
    for a=1:size(LRR,1)
        for b=1:size(LRR,2)
            ixs_corrected = LRR{a,b}<minrate;
            LRR{a,b}(ixs_corrected) = minrate;
            LRRcorrected(a,b) = sum(ixs_corrected);
        end
    end

    LRRs = arrayfun(@(c){cell2mat(LRR(:,c))},(1:size(LRR))');

    %average/stdev over all LRRs
    LRRall = vertcat(LRRs{:});  

    if isfield(x,'distances')
        G.dists{i} = x.distances;
       
        if ~isnan(zthresh)
            dd = G.dists{i};
            dd(dd==0) = nan;
            avg_dist = nanmean(dd);
            z_avg_dist = zscore(avg_dist);

            G.ix_out{i} = z_avg_dist>zthresh;
        end

    end    

    if ismember('ix_out',G.Properties.VariableNames)
        LRR_ixs = combvec(1:G.numel_n(i),1:G.numel_n(i))';
        ix_out = find(G.ix_out{i});
        pix_out = ismember(LRR_ixs(:,1),ix_out) | ismember(LRR_ixs(:,2),ix_out);
        LRRall = LRRall(~pix_out,:); 
        G.n_out(i) = length(ix_out);

        LRRs = cellfun(@(c){c(~G.ix_out{i},:)},LRRs);
        LRRs = LRRs(~G.ix_out{i});
    end    
    
    %average/stdev over test signals (rows) for each reference signal (col)
    G.lrrs_gmean{i} = cellfun(@(c){geom_mean(c)},LRRs');
    G.lrrs_gsd{i} = cellfun(@(c){geom_std(c)},LRRs');

    G.lrrs_amean{i} = cellfun(@(c){nanmean(c)},LRRs'); %#ok<*NANMEAN> 
    G.lrrs_asd{i} = cellfun(@(c){nanstd(c)},LRRs');    %#ok<*NANSTD> 

    G.lrr_gmean(i,:) = geom_mean(LRRall);
    G.lrr_gsd(i,:) = geom_std(LRRall);    

    G.lrr_amean(i,:) = nanmean(LRRall);
    G.lrr_asd(i,:) = nanstd(LRRall); 

    G.lrr_corrected{i} = LRRcorrected;

end
status('reset');

%%

G.onset = cellfun(@(c){c(1)},G.cond);
G.coda = cellfun(@(c){c(end)},G.cond);

subjs = unique(G.subj);
onsets = unique(G.onset);
codas = unique(G.coda);

G.subjix = cellfun(@(c)find(ismember(subjs,c)),G.subj);
G.onsetix = cellfun(@(c)find(ismember(onsets,c)),G.onset);
G.codaix = cellfun(@(c)find(ismember(codas,c)),G.coda);

G.t = repmat({Gi.t},height(G),1);

%% summarize exclusions

G.n_corrected = cellfun(@(c)sum(c>0,'all'),G.lrr_corrected);
fprintf('\t%d LRRs corrected for low/negative slope\n',sum(G.n_corrected));
fprintf('\t%d/%d (%1.3f%%) tokens excluded as outliers\n',sum(G.n_out),sum(G.numel_n),sum(G.n_out)/sum(G.numel_n));

%%

save(outfile,'G','Gi');


end
