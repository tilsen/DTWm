function [map,Xali] = ltwm(X,varargin)
%ltwm applies linear time warping to a set of 1-d signals, which are input in a
%cell array, e.g.
%
%  [map,Xali] = ltwm(X)
%
%  [map,Xali] = ltwm(X,'length',lengthopt)
%               lengthopt is a target length method: median
%               (default), mean, mode, or an integer

p = inputParser;

default_length = 'median';

validinput = @(c)(isnumeric(c) & numel(c)>1) | iscell(X);

addRequired(p,'X',validinput);
addParameter(p,'length',default_length);

parse(p,X,varargin{:});

r = p.Results;

if ischar(r.length)
    switch(r.length)
        case 'median'
            L = median(cellfun(@(c)length(c),X));
        case 'mean'
            L = round(mean(cellfun(@(c)length(c),X)));
        case 'mode'
            L = mode(cellfun(@(c)length(c),X));
        otherwise
            L = median(cellfun(@(c)length(c),X));
    end
elseif isnumeric(r.length)
    L = round(r.length);
end

if ~iscell(X)
    X = {X};
end

for i=1:numel(X)
    
    x = X{i};
    Lx = length(x);
    io = 1:Lx;
    ii = linspace(1,Lx,L);

    Xali{i} = interp1(io,x,ii);

    map{i} = ii;

end



end

