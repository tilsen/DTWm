function [varargout] = stf(varargin)

dbstop if error;
p = inputParser;

if ~exist('inmarg','var'), inmarg=[]; end
if ~exist('exmarg','var'), exmarg=[]; end
if ~exist('axpan','var'), axpan=[]; end

default_axpan = [1 1];
default_aspect = Inf;
default_bare = true;
default_exmarg = [0.05 0.05 0.01 0.01];
default_inmarg = [0.05 0.05];
default_handlearray = 'vector';

addOptional(p,'axpan',default_axpan,@(x)isnumeric(x) || ishandle(x));
addOptional(p,'exmarg',default_exmarg,@(x)isnumeric(x));
addOptional(p,'inmarg',default_inmarg,@(x)isnumeric(x));
addParameter(p,'aspect',default_aspect);
addParameter(p,'bare',default_bare);
addParameter(p,'handlearray',default_handlearray);

%parse(p,axpan,exmarg,inmarg,varargin{:});
parse(p,varargin{:});

make_bare = @(h)set(h,'menubar','none','toolbar','none','name','','numbertitle','off');

defpos = get(0,'defaultFigurePosition');

if nargin==0
    figh = figure('units','normalized','position',defpos);
    if p.Results.bare, make_bare(figh); end
    setaspect(p.Results.aspect,figh);
    varargout = {figh};
    return;
end

figh = figure('units','normalized','position',defpos);
if p.Results.bare, make_bare(figh); end

if isinf(p.Results.aspect)
    set(figh,'WindowState','maximized');
else
    setaspect(p.Results.aspect,figh);
end

ax = stfig_axpos(p.Results.axpan,[p.Results.exmarg p.Results.inmarg]);

switch(p.Results.handlearray)
    case 'matrix'
        if all(size(p.Results.axpan)==[1 2])
            ax = reshape(ax,p.Results.axpan(2),[])';
        else
            ax = ax(p.Results.axpan);
        end
end

varargout = {ax,figh};

end

function [axh,varargout] = stfig_axpos(axn,margins)
% generate axes, using normalized figure coordinates


if nargin<2, margins = [.08 .08 .01 .02 .005 .005]; end %default margins
if length(margins)~=6
    if length(margins)==4
        margins(5:6) = [0.005 0.005];
    else
        fprintf('error: must specify 4 or 6 margin values\n'); return;
    end
end

imarg = margins(5:6);
xmarg = margins([1 3]);
ymarg = margins([2 4]);

%%
if(all(size(axn)==[1 2]))  %convert [rows, cols] to grid specification
    
    nr = axn(1);
    nc = axn(2);
    axn = reshape((1:(nr*nc)),nc,nr)';   
end

%%

axn = flipud(axn);
nr = size(axn,1);
nc = size(axn,2);
arng = 1 - sum(xmarg) + imarg(1);
brng = 1 - sum(ymarg) + imarg(2);

%convert to offsets
w = arng/nc;
h = brng/nr;

xoff = xmarg(1):w:(1-xmarg(2));
yoff = ymarg(1):h:(1-ymarg(2));

if length(xoff)<size(axn,2)
   xoff = linspace(xmarg(1),(1-xmarg(2)),size(axn,2)); w = mean(diff(xoff));
end
if length(yoff)<size(axn,1)
   yoff = linspace(ymarg(1),(1-ymarg(2)),size(axn,1)); h = mean(diff(yoff));
end
%---------

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

set(axh,'ActivePositionProperty','position');

end

function [varargout] = setaspect(varargin)

dbstop if error;
p = inputParser;

default_figh = nan;
default_monitornum = 1;

monpos = get(0,'MonitorPositions');
maxasp = @()monpos(1,3)/monpos(1,4);
default_aspect = inf;

addOptional(p,'target_aspect',default_aspect,@(x)isnumeric(x));
addOptional(p,'figh',default_figh,@(x)ishandle(x));
addOptional(p,'monitornum',default_monitornum,@(x)isnumeric(x));
varargout{1} = [];

parse(p,varargin{:});

figh = p.Results.figh;
if ~ishandle(figh)
    figh = gcf;
end

target_aspect = p.Results.target_aspect;
if isinf(target_aspect)
    %target_aspect =  maxasp();
    set(gcf,'WindowState','maximized');
    return;
end

set(figh,'units','inches');
set(gcf,'WindowState','maximized'); drawnow;

newpos = figh.OuterPosition;
curr_aspect = newpos(3)/newpos(4);

adjfac = target_aspect/curr_aspect;

if curr_aspect > target_aspect %too wide, reduce width
    newpos(3) = newpos(3)*adjfac;
elseif curr_aspect < target_aspect %reduce height
    newpos(4) = newpos(4)/adjfac;
end

set(figh,'outerposition',newpos); drawnow;
set(figh,'units','normalized');

end
