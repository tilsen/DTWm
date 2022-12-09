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
def_colormap = parula(100);
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

%%
function [axh] = stfig_subaxpos(ax,axn,margins)

%ax: axes handle or array of axes handles
%axn: new layout spec or cell array of new layout specs
%margins: new margins (relative to parent) or cell array of new margins

if numel(ax)>1
    ax = ax(:);
    if ~iscell(axn)
        axn = repmat({axn},1,length(ax));
    end
    if ~iscell(margins)
        margins = repmat({margins},1,length(ax));
    end    
    
    for i=1:length(ax)
        axh{i,1} = stfig_subaxpos(ax(i),axn{i},margins{i});
    end
    return
end

pos = ax.Position;
pos_l = pos(1);
pos_r = pos(1)+pos(3);
pos_b = pos(2);
pos_t = pos(2)+pos(4);

if all(size(axn)==[1 2]) %convert to grid
    nr = axn(1);
    nc = axn(2);
    axn = reshape((1:(nr*nc)),nc,nr)';
end

if nargin<3, margins = [.01 .01 .01 .01 .01 .01]; end %default margins

imarg = margins(5:6);
xmarg = margins([1 3]);
ymarg = margins([2 4]);

axn = flipud(axn);
nr = size(axn,1);
nc = size(axn,2);

xoff = linspace(pos_l+xmarg(1),(pos_r-xmarg(2)),size(axn,2)+1); w = mean(diff(xoff));
yoff = linspace(pos_b+ymarg(1),(pos_t-ymarg(2)),size(axn,1)+1); h = mean(diff(yoff));

%----------
nax = length(unique(axn(~isnan(axn(:)))));
axposc = cell(1,nax);
for y=1:nr
    for x=1:nc
        if isnan(axn(y,x)), continue; end
        axposc{axn(y,x)} = [axposc{axn(y,x)}; xoff(x) yoff(y) xoff(x)+(w-imarg(1)) yoff(y)+(h-imarg(2))];
    end
end

for j=1:length(axposc)
    axpos(j,:) = [min(axposc{j}(:,1)) min(axposc{j}(:,2)) max(axposc{j}(:,3)) max(axposc{j}(:,4))]; %#ok<AGROW>
    
end
axpos(:,3:4) = axpos(:,3:4)-axpos(:,1:2);

varargout{1} = axpos;
for i=1:size(axpos,1)
    axh(i) = axes('position',axpos(i,:)); %#ok<AGROW>
end

end

function [h] = stfig_gridlines(xy,varargin)

p = inputParser;

def_xy = 'xy';
def_parent = gca;
def_color = [0 0 0];
def_linewidth = 1;
def_xvals = nan;
def_yvals = nan;

addOptional(p,'xy',def_xy,@(x)ischar(x));
addParameter(p,'parent',def_parent);
addParameter(p,'color',def_color);
addParameter(p,'linewidth',def_linewidth);
addParameter(p,'xvals',def_xvals);
addParameter(p,'yvals',def_yvals);

parse(p,xy,varargin{:});

axh = p.Results.parent;

X = p.Results.xvals(:)';
if isnan(X)
    X = min(axh.XLim):max(axh.XLim);  
    Xr = axh.XLim;
else
    Xr = minmax(X);
end

nX = length(X);

Y = p.Results.yvals(:)';
if isnan(Y)
    Y = min(axh.YLim):max(axh.YLim);
    Yr = axh.YLim;
else
    Yr = minmax(Y);
end
nY = length(Y);

h = {};
if contains(p.Results.xy,'x')
   h{end+1} = line( repmat(X,2,1), repmat(Yr',1,nX),...
       'color',p.Results.color,'parent',axh,'linew',p.Results.linewidth);
end
if contains(p.Results.xy,'y')
   h{end+1} = line( repmat(Xr',1,nY), repmat(Y,2,1),...
       'color',p.Results.color,'parent',axh,'linew',p.Results.linewidth);
end

end

function [th] = matrix_text(M,varargin)

dbstop if error;

p = inputParser;

def_fontsize = get(0,'defaultaxesFontSize');
def_fontcolor = get(0,'defaultaxesColor');
def_fontweight = 'normal';
def_color_threshold = nan;
def_formatstr = '%1.2f';
def_parent = nan;
def_x = nan;
def_y = nan;

addRequired(p,'M',@(x)ismatrix(x) | iscell(x));
addParameter(p,'fontsize',def_fontsize);
addParameter(p,'fontcolor',def_fontcolor);
addParameter(p,'fontweight',def_fontweight);
addParameter(p,'color_threshold',def_color_threshold);
addParameter(p,'formatstr',def_formatstr);
addParameter(p,'parent',def_parent);
addParameter(p,'x',def_x);
addParameter(p,'y',def_y);

parse(p,M,varargin{:});

res = p.Results;

if ~ishandle(res.parent)
    res.parent = gca;
end

if ismatrix(M)
    Mstr = arrayfun(@(c)sprintf(res.formatstr,c),M,'un',0);
    Mstr(isnan(M)) = {''};
else
    Mstr = M;
end

if ~isnan(res.color_threshold) && size(res.fontcolor,1)==1
    res.fontcolor = repmat(res.fontcolor,2,1);
end

if isnan(res.x)
    res.x = 1:size(M,2);
end
if isnan(res.y)
    res.y = 1:size(M,1);
end

for r=1:size(M,1)
    for c=1:size(M,2)

        
        color = res.fontcolor(1,:);  
        
        if ~isnan(res.color_threshold)
            if M(r,c)<=res.color_threshold
                color = res.fontcolor(1,:);
            else
                color = res.fontcolor(2,:);
            end       
        end

        th(r,c) = text(res.x(c),res.y(r),Mstr{r,c},...
            'hori','center', ...
            'verti','mid', ...
            'color',color, ...
            'fontsize',res.fontsize, ...
            'fontweight',res.fontweight, ...
            'parent',res.parent);
    end
end


end


function [] = text_color_threshold(imh,th,varargin)

p = inputParser;

def_colors = [0 0 0; 1 1 1];
def_thresholds = 127.5;

addRequired(p,'imh',@(x)ishandle(imh) && numel(imh.CData)==numel(th));
addRequired(p,'th',@(x)all(ishandle(th),'all'));
addParameter(p,'thresholds',def_thresholds);
addParameter(p,'colors',def_colors);

parse(p,imh,th,varargin{:});

r = p.Results;

if size(r.colors,1)~=(numel(r.thresholds)+1)
    fprintf('error: number of colors must equal number of thresholds + 1\n');
    return;
end

cmap = get(imh.Parent,'Colormap');
cdata = imh.CData;
clim = get(imh.Parent,'CLim');

for i=1:numel(cdata)
    sc = (cdata(i)-clim(1))/range(clim);
    ix = min(1+floor(sc*size(cmap,1)),size(cmap,1));
    colors(i,:) = cmap(ix,:);
end

colors = colors*255;

f = [.299 .587 .114];
hsp = sqrt(sum(f.*colors.^2,2));

set(th,'color',r.colors(1,:));
for i=1:length(r.thresholds)
    set(th(hsp>r.thresholds(i)),'Color',r.colors(i+1,:));
end

end


