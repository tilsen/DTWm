function [] = lrr_inspect(X,Y,varargin)
%function for inspecting lrr and dtwm together
%
%inputs: X, Y - signals to align
%        

p = inputParser;

addRequired(p,'X');
addRequired(p,'Y');
addParameter(p,'dtw_params',{});
addParameter(p,'lrr_params',{});

parse(p,X,Y,varargin{:});

res = p.Results;

[map,distxy,info] = dtwm(X,Y, res.dtw_params{:});

%slopes/LRRs of columns (targets) to rows (comparisons)
[local_slope,singular_warning] = lrr(map, res.lrr_params{:});

h=[];
h = setup_gui(h,res,info,local_slope);

% while 1
%     uiwait;
% 
% end


end

%%
function [h] = setup_gui(h,res,info,local_slope)

if isempty(h)
    h.ax = stf([1 1 2; 1 1 nan],[0.05 0.05 0.01 0.25],[0.10 0]);
end

step_pattern = res.dtw_params{find(ismember(res.dtw_params,'step_pattern'))+1};


h.dtwm = plot_dtw_matrix(info,res.X,res.Y,'parent',h.ax(1));
drawnow;
adjwidth(h.dtwm.cbh,-0.025);
shiftposx(h.dtwm.cbh,-0.015);

h.steppat = plot_step_pattern(select_step_pattern(step_pattern),'parent',h.ax(2));

pos = h.dtwm.axd.Position + [0 0.01 0 0];

h.ax(3) = axes('Position',[pos(1) sum(pos([2 4])) pos(3) 0.99-sum(pos([2 4]))]);
plot(local_slope,'k-','linew',2);
axis tight;
axrescaley(0.05);

set(h.ax(3),'XGrid','on','YGrid','on','tickdir','out','Box','off');

end