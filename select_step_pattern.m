function [s] = select_step_pattern(pattern)

s = def_step_patterns;

if nargin>0 
    if ~iscell(pattern)
        pattern = {pattern};
    end

    ix_not_specified = find(~ismember(pattern,s.pattern));
    if ~isempty(ix_not_specified)
        fprintf('The following step patterns are not specified:\n');
        arrayfun(@(c)fprintf('%s\n',pattern{c},ix_not_specified));
        pattern = pattern(setdiff(1:length(pattern),ix_not_specified));
    end

    ixs = cellfun(@(c)find(ismember(s.pattern,c)),pattern);
    s = s(ixs,:);
end

if height(s)==1
    s = table2struct(s);
end


end