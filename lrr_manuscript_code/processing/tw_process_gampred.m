function [] = tw_process_gampred(model)

dbstop if error;
h = tw_helpers;

warning('off','MATLAB:table:ModifiedAndSavedVarnames');

if nargin==0
    files = rdir([h.gam_dir 'gampred_*.csv']);
else
    files = rdir([h.gam_dir 'gampred_' model '.csv']);
end

%% to reconstruct condition index numbers:
segs = {'0','p','t'};
conds = strcat(repmat(segs,numel(segs),1),'_',repmat(segs',1,numel(segs)));
conds = conds(:);

%%
for i=1:length(files)
    [~,f] = fileparts(files(i).name);

    status_str = status(sprintf('%s: %i/%i',f,i,length(files))); %#ok<NASGU> 

    x = readtable(files(i).name,'TreatAsMissing',{'NA'});

    D = x(~isnan(x.cond),:);
    D.Var1 = [];
    D.subj = [];
    D = sortrows(D,{'cond' 'time'});

    ix_all_nan = all(isnan(table2array(D)));
    D = D(:,~ix_all_nan);

    E = x(~isnan(x.d_1_2),:);
    E.Var1 = [];
    ix_all_nan = all(isnan(table2array(E)));
    E = E(:,~ix_all_nan);

    E = sortrows(E,{'time'});

    [~,f] = fileparts(files(i).name);

    Ds = D;
    Dd = E;

    Ds = fits_by_rows(Ds,conds);
    Dd = fits_by_rows(Dd,conds);

    Ds.name = repmat(f,height(Ds),1);
    Dd.name = repmat(f,height(Dd),1);

    if contains(f,'log')

        %b/c get_pred returns 1.96*s.e. in log space here for the fixed effect:
        ixa = ismember(Ds.subj,'ALL');
        Ds.ci(ixa,:) = (exp(Ds.mu(ixa,:)+Ds.ci(ixa,:)/norminv(0.975))-exp(Ds.mu(ixa,:))); 
        Ds.ci(~ixa,:) = (exp(Ds.mu(~ixa,:)+Ds.ci(~ixa,:))-exp(Ds.mu(~ixa,:))); 
        Ds.mu = exp(Ds.mu);
        
        ixa = ismember(Dd.subj,'ALL');
        Dd.ci(ixa,:) = (exp(Dd.mu(ixa,:)+Dd.ci(ixa,:)/norminv(0.975))-exp(Dd.mu(ixa,:))); 
        Dd.ci(~ixa,:) = (exp(Dd.mu(~ixa,:)+Dd.ci(~ixa,:))-exp(Dd.mu(~ixa,:))); 

    else
        %convert to 1.0 s.e.
        ixa = ismember(Ds.subj,'ALL');
        Ds.ci(ixa,:) = Ds.ci(ixa,:)/norminv(0.975);

        ixa = ismember(Dd.subj,'ALL');
        Dd.ci(ixa,:) = Dd.ci(ixa,:)/norminv(0.975);        

    end

    save([h.gam_dir f '.mat'],'Ds','Dd');

end
status('reset');

end

%%
function [D] = fits_by_rows(X,conditions)

vn = X.Properties.VariableNames;

switch(ismember('cond',vn))
    case 1
        conds = unique(X.cond);
        traj = vn(contains(vn,'mu'));
        D=[];
        for i=1:length(conds)
            ix = ismember(X.cond,conds(i));
            for j=1:length(traj)
                D(end+1).name = traj{j};
                D(end).cond = conds(i);
                D(end).mu = X.(traj{j})(ix)';
                D(end).ci = X.(strrep(traj{j},'mu','ci'))(ix)';
                D(end).t = X.time(ix)';
            end
        end
        D = struct2table(D);
        D.cond = conditions(D.cond);
        D.subj = strrep(D.name,'mu_','');

        subjn = cellfun(@(c)str2double(c),regexp(D.name,'\d+$','match','once'));
        D.subj(~isnan(subjn)) = arrayfun(@(c){sprintf('P%02d',c)},subjn(~isnan(subjn)));
        D.subj(strcmp(D.subj,'mu')) = {'ALL'};

    case 0

        %subject-indepedent pred:
        traj = regexp(vn,'^d_\d+_\d+$','match','once');
        traj = traj(~cellfun('isempty',traj));

        D=[];
        for i=1:length(traj)
            D(end+1).name = traj{i};
            D(end).subj = 'ALL';
            D(end).cond = strrep(traj{i},'d_','');
            D(end).mu = X.(traj{i})';
            D(end).ci = X.(strrep(traj{i},'d_','ci_'))';
            D(end).t = X.time';
        end

        %subject-specific pred:
        traj = regexp(vn,'^d_\d+_\d+_\d+$','match','once');
        traj = traj(~cellfun('isempty',traj));
        for i=1:length(traj)
            D(end+1).name = traj{i};
            parts = regexp(traj{i},'_','split');
            D(end).subj = parts{end};
            D(end).cond = [parts{2} '_' parts{3}];
            D(end).mu = X.(traj{i})';
            D(end).ci = X.(strrep(traj{i},'d_','ci_'))';
            D(end).t = X.time';
        end        

        D = struct2table(D);

        for i=1:height(D)
            if ~contains(D.subj{i},'ALL')
                D.subj{i} = sprintf('P%02d',str2double(D.subj{i}));
            end
        end

        Dr = D;
        Dr.mu = -Dr.mu;

        Dr.cond = cellfun(@(c){fliplr(c)},Dr.cond);
        
        D = [D; Dr];

        D.cond1 = conditions(cellfun(@(c)str2double(c(1)),D.cond));
        D.cond2 = conditions(cellfun(@(c)str2double(c(end)),D.cond));
end

end
