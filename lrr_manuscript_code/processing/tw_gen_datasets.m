function [] = tw_gen_datasets()

dbstop if error;

h = tw_helpers;

overwrite = true;

%% 
%MFCC parameters
mfccpars.n_coeffs = 13;
mfccpars.winsize = 0.025;
mfccpars.stepsize = 0.001;

%anchorpoint and time horizon
anchf = 'V_vel';
tlims = [-0.750 0.750];

%datasets:

S(1).dataset = 'CVCarl';
S(end).signals = 'articulatory';
S(end).trng = tlims;

S(end+1).dataset = 'CVCacl';
S(end).signals = 'acoustic';
S(end).trng = tlims;

S = struct2table(S);

%%
X = load([h.data_dir 'onsetcoda_opendata.mat']);
X.D = recenter_trim_articulatory_signals(X.D,anchf,tlims);
X.D = recenter_trim_acoustic_signals(X.D,anchf,tlims);

X.D = add_mfcc(X.D,mfccpars);
X.D.wav = [];

%%
for i=1:height(S)

    outfile = [h.datasets_dir 'dataset_' S.dataset{i} '.mat'];

    if ~overwrite && exist(outfile,'file'), continue; end

    %makes data tables
    %UserData{1}: table of groups
    %UserData{2}: data table indexing function (by group id)
    %UserData{3}: additional dataset-specific metadata; 

    fprintf('generating dataset: %s\n',S.dataset{i});

    D = collect_data(S(i,:),X.D);
    D = limit_times(D,S.trng(i,:));

    save(outfile,'D');
end


end

%%
function [D] = collect_data(S,D)

D.onset = cellfun(@(c){c(1)},D.cond);
D.coda = cellfun(@(c){c(end)},D.cond);

metadata = D.Properties.UserData{1};

switch(S.signals{1})
    case 'articulatory'
        D.X = D.kin;          
        metadata.t = D.Properties.UserData{1}.t_kin;
        metadata.chans = metadata.kin_chans;
        D.kin = []; D.mfcc = [];
        
    case 'acoustic'
        D.X = D.mfcc;
        metadata.t = metadata.t_mfcc;
        D.kin = []; D.mfcc = [];
        
end

if contains(S.dataset,'CVX')
    D.cond = D.onset;
end

if contains(S.dataset,'XVC')
    D.cond = D.coda;
end

ix_out = D.resp_error | D.wav_outofrange | D.kin_outofrange;
fprintf('tokens excluded: %d\n',sum(ix_out));

keepf = {'subj' 'fname' 'cond' 'onset' 'coda' 'X'};
D = D(~ix_out,keepf);

D = sortrows(D,{'subj','cond'});

%add metadata
D.n = true(height(D),1);
G = grpstats(D,{'subj' 'cond'},'numel','DataVars','n');
G.id = (1:height(G))';
D.n = [];

D.Properties.UserData{1} = G;
D.Properties.UserData{2} = @(gid)ismember(D.subj,G.subj(gid)) & ismember(D.cond,G.cond(gid));
D.Properties.UserData{3} = metadata;

end

%% recenter articulatory signals
function [D] = recenter_trim_articulatory_signals(D,anchf,tlims)

info = D.Properties.UserData{1};
kinFs = info.kin_Fs;

D.t_kin = cellfun(@(c,d){c-d},D.t_kin,num2cell(D.(anchf)));

% mark out of range time coordinates
D.kin_outofrange = cellfun(@(c)c(1)>tlims(1) | c(end)<tlims(2),D.t_kin);

t_kin = tlims(1):(1/kinFs):tlims(2);
len = length(t_kin);

X = nan(height(D),length(info.kin_chans),len);
for i=1:height(D)
    ix = find(D.t_kin{i}>=tlims(1),1,'first');
    X(i,:,:) = D.kin{i}(:,ix+(0:len-1));
end

D.kin = single(X);

D.t_kin = [];
D.Properties.UserData{1}.t_kin = t_kin;

end

%% recenter acoustic signals
function [D] = recenter_trim_acoustic_signals(D,anchf,tlims)

len_wav = size(D.wav,2);
wavFs = D.Properties.UserData{1}.wav_Fs;
t_wav = (0:(len_wav-1))/wavFs;

D.t_wav = repmat(t_wav,height(D),1);
D.t_wav = D.t_wav-D.(anchf);

%desired plus extra to accomodate mfcc window:
tlims = tlims + [-0.025 0.025]; 

% mark out of range time coordinates
D.wav_outofrange = D.t_wav(:,1)>tlims(1) | D.t_wav(:,end)<tlims(2);

len = round(diff(tlims)*wavFs);

X = nan(height(D),len);
for i=1:height(D)
    ix = find(D.t_wav(i,:)>=tlims(1),1,'first');
    X(i,:) = D.wav(i,ix:ix+len-1);
end

t_wav = tlims(1):(1/wavFs):tlims(2);

D.wav = X;
D.t_wav = [];
D.Properties.UserData{1}.t_wav = t_wav;

end

%% add mfccs
function [D] = add_mfcc(D,mfccpars)

audio_Fs = D.Properties.UserData{1}.wav_Fs;

mfccpars.winlen = round(audio_Fs*mfccpars.winsize);
mfccpars.n_overlap = round((mfccpars.winsize-mfccpars.stepsize)*audio_Fs);
mfccpars.window = hamming(mfccpars.winlen,'periodic');

D.Properties.UserData{1}.mfcc_chans = arrayfun(@(c){sprintf('mfcc_%02d',c)},(1:mfccpars.n_coeffs)');
D.Properties.UserData{1}.mfcc_Fs = 1/mfccpars.stepsize;
D.Properties.UserData{1}.mfcc_dt = mfccpars.stepsize;

for i=1:height(D)

    status_str = status(sprintf('processing mfcc %05i/%i',i,height(D))); %#ok<NASGU> 

	x = D.wav(i,:)';
    
	[coeffs,~,~,locs] = mfcc(x,audio_Fs,...
        'Window',mfccpars.window, ...
        'OverlapLength',mfccpars.n_overlap, ...
        'NumCoeffs',mfccpars.n_coeffs,...
        'LogEnergy','ignore');

    D.mfcc{i} = single(coeffs');

end
status('reset');

D.mfcc = permute(cat(3,D.mfcc{:}),[3 1 2]);

t_mfcc = locs'/audio_Fs - mfccpars.winsize/2 + D.Properties.UserData{1}.t_wav(1);

D.Properties.UserData{1}.t_mfcc = t_mfcc;

end

%%
function [D] = limit_times(D,trng)
t = D.Properties.UserData{end}.t;
ixs = t>=trng(1) & t<=trng(2);
D.Properties.UserData{end}.t = t(ixs);
D.X = D.X(:,:,ixs);
end

