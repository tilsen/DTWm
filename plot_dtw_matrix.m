function [h] = plot_dtw_matrix(dtwinfo,s1,s2,varargin)
% plot_dtw_matrix
%
%   plots a distance matrix from dtwm with signals in margins 

dbstop if error;
p = inputParser;

def_parent = nan;
def_path_linespec = '-';
def_path_linecolor = 'w';
def_plot_signals = true;
def_signals_linespec = {'-','-'};
def_colormap = viridis(100);
def_gridsize = 20;
def_matrixsize = 16;
def_scolors = lines(2);
def_matrixgrid = false;
def_matrixfield = 'localCostMatrix';
def_showvalues = false;
def_formatstr = '%1.2f';
def_fontsize = get(0,'defaulttextfontSize');
def_fontcolor = [0 0 0];
def_fontweight = 'normal';

addRequired(p,'dtwinfo',@(x)isstruct(x));
addOptional(p,'s1',@(x)isvector(x));
addOptional(p,'s2',@(x)isvector(x));
addParameter(p,'parent',def_parent);
addParameter(p,'matrixgrid',def_matrixgrid);
addParameter(p,'path_linespec',def_path_linespec);
addParameter(p,'path_linecolor',def_path_linecolor);
addParameter(p,'plot_signals',def_plot_signals);
addParameter(p,'signals_linespec',def_signals_linespec);
addParameter(p,'colormap',def_colormap);
addParameter(p,'gridsize',def_gridsize);
addParameter(p,'matrixsize',def_matrixsize); %width of signal axes in grid units
addParameter(p,'scolors',def_scolors);
addParameter(p,'matrixfield',def_matrixfield);
addParameter(p,'showvalues',def_showvalues);
addParameter(p,'formatstr',def_formatstr);
addParameter(p,'fontsize',def_fontsize);
addParameter(p,'fontcolor',def_fontcolor);
addParameter(p,'fontweight',def_fontweight);

parse(p,dtwinfo,s1,s2,varargin{:});

r = p.Results;

if ~ishandle(r.parent)
    ax = gca;
else
    ax = r.parent;
end

if r.plot_signals
    axg = nan(r.gridsize);
    ms = r.matrixsize;
    axg(1:ms,(end-ms+1):end) = 3;
    axg(1:ms,1:(end-ms-1)) = 2;
    axg(ms+2:end,(end-ms+1):end) = 1;
    axsub = stfig_subaxpos(ax,axg,[0 0 0 0 0 0]);
    axd = axsub(3);
    axs = axsub(1:2);
    delete(ax);
else
    axd = ax;
    axs = [];
end

if ~isfield(dtwinfo,r.matrixfield)
    fprintf('error: matrix field not found\n'); return;
end

M = dtwinfo.(r.matrixfield)';

%plot local cost matrix and warping path
h.imh = imagesc(M,'parent',axd); 
hold(axd,'on');
colormap(axd,r.colormap);
set(axd,'YDir','normal');
h.cbh = colorbar(axd);
h.cbh.Position(1) = sum(axd.Position([1 3]))+0.025;
h.cbh.Position(3) = 0.05;

h.maph = plot(dtwinfo.index1s,dtwinfo.index2s,r.path_linespec,...
    'color',r.path_linecolor,...
    'parent',axd);

if r.matrixgrid
    h.gridh = stfig_gridlines('xy');
end

if r.plot_signals
    set(axd,'YTickLabel',[],'XTickLabel',[]);
    h.s1 = plot(r.s1',...
        r.signals_linespec{1},...
        'parent',axs(1),...
        'color',r.scolors(1,:));
    h.s2 = plot(r.s2,(1:length(r.s2)),...
        r.signals_linespec{2},...
        'parent',axs(2), ...
        'color',r.scolors(end,:));
    hold(axs,'on');
    axis(axs,'tight');
    ylim(axs(1),ylim(axs(1))+0.05*diff(ylim(axs(1)))*[-1 1]);
    xlim(axs(1),axd.XLim);
    xlim(axs(2),xlim(axs(2))+0.05*diff(xlim(axs(2)))*[-1 1]);
    ylim(axs(2),axd.YLim);
    set(axs,'XGRid','on','ygrid','on');
    h.xlabs1 = xlabel(axs(1),'sample index');
    h.ylabs2 = ylabel(axs(2),'sample index');
    set(axs(2),'XDir','reverse');
    axs(1).XTick = unique([1 axs(1).XTick]);
    axs(2).YTick = unique([1 axs(2).YTick]);
end

if r.showvalues
    h.vh = matrix_text(M,'formatstr',r.formatstr,'fontsize', ...
        r.fontsize,'fontcolor',r.fontcolor, ...
        'fontweight',r.fontweight);

    text_color_threshold(h.imh,h.vh, ...
        'colors',[1 1 1; 0 0 0]);
end

h.axd = axd;
h.axs = axs;

end

