function [h] = plot_step_pattern(pat,varargin)

p = inputParser;

%defaultgrid = false;
def_linewidth = 2;
def_fontsize = 16;
def_parent = nan;
def_markersize = 4;

addRequired(p,'pat',@(x)isstruct(x) | ismatrix(x));
addParameter(p,'linewidth',def_linewidth);
addParameter(p,'parent',def_parent);
addParameter(p,'fontsize',def_fontsize);
addParameter(p,'markersize',def_markersize);

parse(p,pat,varargin{:});

r = p.Results;

if isstruct(pat)
    M = pat.mx;
    titlestr = pat.pattern;
    if ~isempty(pat.norm)
        titlestr = [titlestr ' (' pat.norm{:} ')'];
    end
else
    M = pat;
    titlestr = '';
end

if isempty(r.parent) || ~ishandle(r.parent)
    r.parent = gca;
end

h = struct('ph',[],'th',[]);

ids = unique(M(:,1));

%rightward/upward orientation
M(:,2:3) = -M(:,2:3);

for i=1:length(ids)
    ix = M(:,1)==ids(i);
    m = M(ix,:);

    xx = m(:,3);
    yy = m(:,2);

    costs = m(2:end,4);

    h.ph(end+1) = plot(xx,yy,'ko-',...
        'parent',r.parent,'markerfacecolor','w',...
        'linewidth',r.linewidth,'markersize',r.markersize); 

    hold(r.parent,'on');
    
    plot(xx(1),yy(1),'ko','markerfacecolor','k','markersize',r.markersize);

    for j=1:length(costs)
        
        if mod(costs(j),1)==0
            frmtstr = '%1.0f';
        else
            frmtstr = '%1.2f';
        end

        xv = xx([j j+1]);
        yv = yy([j j+1]);

        astrs = alignments(xv,yv);

        h.th(end+1) = text(mean(xv),mean(yv),...
            num2str(costs(j),frmtstr),...
            astrs{:},'parent',r.parent,...
            'fontsize',r.fontsize,'margin',0.001);
    end
end

h.phc = plot(0,0,'ks','Markerfacecolor','k',...
    'parent',r.parent,'markersize',r.markersize+2);

h.titleh = title(r.parent,titlestr);

axis(r.parent,'tight');
lims = minmax([xlim(r.parent) ylim(r.parent)]);
N = min(lims);

set(r.parent,'XTick',[N:0],'YTick',[N:0]);

lims = lims+0.05*diff(lims)*[-1 1];
set(r.parent,'xlim',lims,'ylim',lims,...
    'box','off','tickdir','out','xgrid','on','ygrid','on');


end

%%
function [astrs] = alignments(xv,yv)

astrs = {'hori','left','verti','top'};

if diff(xv)<0.001
    astrs{4} = 'mid';
end
if diff(yv)<0.001
    astrs{2} = 'center';
end


end
