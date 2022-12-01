function [maps,distances,map_infos] = dtwm_pairwise(X,varargin)

p = inputParser;

p.KeepUnmatched = 1;

def_use_parallel = false;

addRequired(p,'X',@(x)iscell(x));
addParameter(p,'use_parallel',def_use_parallel);

parse(p,X,varargin{:});

r = p.Results;

maps = cell(length(X));
distances = nan(length(X));
map_infos = cell(length(X));

if r.use_parallel && isempty(gcp("nocreate"))
    fprintf('warning: ''use_parallel'' option selected, but no parallel pool is running\n');
    r.use_parallel = false;
end

switch(r.use_parallel)
    case true
        for a=1:length(X)
            parfor b=1:length(X)
                [maps{a,b},distances(a,b),map_infos{a,b}] = ...
                    dtwm(X{a},X{b}, varargin{:});
            end
        end
    otherwise
        for a=1:length(X)
            for b=1:length(X)
                [maps{a,b},distances(a,b),map_infos{a,b}] = ...
                    dtwm(X{a},X{b}, varargin{:});
            end
        end
end

end