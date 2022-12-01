function [h] = violin_plot(o,X,varargin)

p = inputParser;

def_orientation = 'vertical';
def_scalefactor = 1;
def_ksdensityargs = {};
def_color = [0 0 0];
def_edgecolor = 'none';
def_facealpha = 0.5;
def_dens = nan;
def_pts = nan;

addRequired(p,'o');
addRequired(p,'X');
addParameter(p,'parent',[]);
addParameter(p,'ksdensityargs',def_ksdensityargs);
addParameter(p,'orientation',def_orientation);
addParameter(p,'scalefactor',def_scalefactor);
addParameter(p,'color',def_color);
addParameter(p,'edgecolor',def_edgecolor);
addParameter(p,'facealpha',def_facealpha);
addParameter(p,'dens',def_dens);
addParameter(p,'pts',def_pts);

parse(p,o,X,varargin{:});

r = p.Results;

if isempty(r.parent)
    r.parent = gca;
end

if isnan(r.dens)
    [f,r.pts] = ksdensity(X,r.ksdensityargs{:});
else
    f = r.dens;
end
if isnan(r.pts)
    r.pts = 1:length(f);
end

pts = r.pts;

sc = r.scalefactor;
if sc~=1
    f = f/max(f);
end

switch(r.orientation)
    case 'vertical'
        h(1) = fill(o+sc*[f -fliplr(f)],[pts fliplr(pts)],r.color, ...
            'EdgeColor',r.edgecolor,'FaceAlpha',r.facealpha);
           
    case 'horizontal'

end

end