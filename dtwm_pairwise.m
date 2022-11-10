function [maps,distances,map_infos] = dtwm_pairwise(X,varargin)

maps = cell(length(X));
distances = nan(length(X));
map_infos = cell(length(X));

for a=1:length(X)
    for b=1:length(X)
        [maps{a,b},distances(a,b),map_infos{a,b}] = ...
            dtwm(X{a},X{b}, varargin{:});
    end
end

end