function [map,Xali] = ltwm(X,varargin)

p = inputParser;

default_length = [];

validinput = @(c)(isnumeric(c) & numel(c)>1) | iscell(X);

addRequired(p,'X',validinput);
addParameter(p,'length',default_length);

parse(p,X,varargin{:});

r = p.Results;

if isempty(r.length)
    L = median(cellfun(@(c)length(c),X));
elseif ischar(r.length)
    switch(r.length)
        case 'median'
            L = median(cellfun(@(c)length(c),X));
        case 'mean'
            L = round(mean(cellfun(@(c)length(c),X)));
        case 'mode'
            L = mode(cellfun(@(c)length(c),X));
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

