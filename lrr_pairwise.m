function [local_slopes,distances] = lrr_pairwise(X,varargin)
% lrr_pairwise conducts a dtw on each pair of input signals and then calculates
% the local relative rate. 
% 
% Input: a one-dimensional cell array of signals. Additional arguments are
% passed to dtwm, e.g.
%
% lrr_pairwise(X,'step_pattern','symmetricP05') - will pass 'step_pattern','symmetricP05' to dtwm.  
%
% Output: for N signals, lrr_pairwise outputs an N x N cell array of lrr
% timeseries, along with a matrix of the corresponding warping distances:
%
% [LRR_array,dists] = lrr_pairwise(X);
%
% IMPORTANT: the output array of this function (i.e. an array of lrr timeseries)
% is organized such that ROWS correspond to target (reference) signals 
% and COLUMNS correspond to comparison (query) signals. 
% Thus the lrr timeseries in each ROW of the output array will have the same length.
%
% Note that in the manuscript: Tilsen & Tiede (2023). Looking within events: 
% examining internal temporal structure with local relative rate, 
% the lrr array is depicted such that COLUMNS correspond to target signals 
% and ROWS correspond to comparison signals. 
% 
% Thus the depiction of the lrr array in the manuscript is a
% transposition of the array that is output by this function. 
% Another way to express the differences is as follows:
%
% This function:
% for a=1:N
%   for b=1:N
%       map = dtwm(X{a},X{b});
%       LRR{a,b} = lrr(map);
%   end
% end
%
% LRR array as shown in the manuscript:
% for a=1:N
%   for b=1:N
%       map = dtwm(X{a},X{b});
%       LRR{b,a} = lrr(map);
%   end
% end


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