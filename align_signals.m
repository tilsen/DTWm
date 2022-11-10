function [a1,a2] = align_signals(map,x_ref,x_test,varargin)

p = inputParser;

default_aligntype = 'expansive';

addRequired(p,'map');
addRequired(p,'x_ref');
addRequired(p,'x_test');
addOptional(p,'aligntype',default_aligntype);

parse(p,map,x_test,x_ref,varargin{:});
r = p.Results;

if any(~ismember(map(1,:),1:length(x_ref)))
    fprintf('ERROR: invalid map\n'); return;
end
if any(~ismember(map(2,:),1:length(x_test)))
    fprintf('ERROR: invalid map\n'); return;
end

switch(r.aligntype)
    case 'expansive'
        a1 = x_ref(map(1,:));
        a2 = x_test(map(2,:));

    case 'compressive'
        map_compressed = map(:,~any(diff([map [0;0]],[],2)==0));
        a1 = x_ref(map_compressed(1,:));
        a2 = x_test(map_compressed(2,:));

    case 'reference'
        ix1 = unique(map(1,:));
        a1 = x_ref(ix1);

        ix2 = arrayfun(@(c)map(2,find(map(1,:)==c,1,'last')),ix1);
        a2 = x_test(ix2);

    otherwise
        fprintf('ERROR: unknown alignment type\n'); return;
end

end
