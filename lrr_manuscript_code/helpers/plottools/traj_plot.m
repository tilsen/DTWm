function [fh,ph] = traj_plot(x,ci,t,varargin)

p = inputParser;

def_parent = gca;
def_color = [0 0 0];
def_edgecolor = 'none';
def_facealpha = 0.25;
def_linewidth = 2;
def_linestyle = '-';

addRequired(p,'x');
addRequired(p,'ci');
addOptional(p,'t',@(c)isvector(c) & length(c)==length(x));
addParameter(p,'parent',def_parent);
addParameter(p,'color',def_color);
addParameter(p,'edgecolor',def_edgecolor);
addParameter(p,'facealpha',def_facealpha);
addParameter(p,'linewidth',def_linewidth);
addParameter(p,'linestyle',def_linestyle);


parse(p,x,ci,t,varargin{:});

r = p.Results;

t = r.t(:);
tt = [t; flipud(t)];

x = x(:);


if any(size(ci)==1)
   ci = ci(:);
   ci = x+[-ci ci];
elseif any(size(ci)==2)
    if size(ci,2)~=2
        ci = ci';
    end    
else
    fprintf('unexpected dimensions of confidence interval timeseries\n');
    return;
end

ci = [ci(:,1); flipud(ci(:,2))];

fh = fill(tt,ci,r.color,'EdgeColor',r.edgecolor,'FaceAlpha',r.facealpha,'parent',r.parent); hold(r.parent,'on');
ph = plot(t,x,'color',r.color,'linew',r.linewidth,'linestyle',r.linestyle,'parent',r.parent); 



end