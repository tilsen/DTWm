function [local_slopes,distances] = lrr_pairwise(X,varargin)

p = inputParser;
p.KeepUnmatched = true;

def_winsize = 0.05;
def_useparallel = false;

addRequired(p,'X');
addParameter(p,'winsize',def_winsize);
addParameter(p,'useparallel',def_useparallel);

parse(p,X,varargin{:});

distances = nan(length(X));

dtwm_inputs = {};
if ~isempty(p.Unmatched)
    fn = fieldnames(p.Unmatched);
    c=1;
    for j=1:length(fn)
        dtwm_inputs{c} = fn{j}; 
        c=c+1;
        dtwm_inputs{c} = p.Unmatched.(fn{j});
    end
end

switch(p.Results.useparallel)
    case false

        local_slopes = cell(length(X),length(X));
        for a=1:length(X)
            for b=1:length(X)
                [map,distances(a,b)] = dtwm(X{a},X{b}, dtwm_inputs{:});

                %slopes/LRRs of columns (targets) to rows (comparisons)
                local_slopes{a,b} = lrr(map,'winsize',p.Results.winsize);
            end
        end

    case true

        for a=1:length(X)
            xa = X{a};
            parfor b=1:length(X)
                [map,distances(a,b)] = dtwm(xa,X{b}, dtwm_inputs{:});
                
                %slopes/LRRs of columns (targets) to rows (comparisons)
                local_slopes{a,b} = lrr(map,'winsize',p.Results.winsize);
            end
        end
end

end