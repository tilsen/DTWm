function [fh,ph] = dens_plot(dens,pnts,varargin)

p = inputParser;

def_parent = gca;
def_color = [0 0 0];
def_edgecolor = 'none';
def_facealpha = 0.25;
def_linewidth = 2;
def_linestyle = '-';
def_orientation = 'horizontal';
def_baseline = 0;
def_outline = 0;

addRequired(p,'dens');
addRequired(p,'pnts');
addOptional(p,'baseline',def_baseline,@(x)isscalar(x));
addParameter(p,'parent',def_parent);
addParameter(p,'color',def_color);
addParameter(p,'edgecolor',def_edgecolor);
addParameter(p,'facealpha',def_facealpha);
addParameter(p,'linewidth',def_linewidth);
addParameter(p,'linestyle',def_linestyle);
addParameter(p,'orientation',def_orientation);
addParameter(p,'outline',def_outline);

parse(p,dens,pnts,varargin{:});

r = p.Results;

pnts = pnts(:);
dens = dens(:);

fr = [dens; r.baseline*ones(size(dens))];
pp = [pnts; flipud(pnts)];

switch(r.orientation)
    case 'horizontal'
        fh = fill(pp,fr,r.color, ...
            'EdgeColor',r.edgecolor,'FaceAlpha',r.facealpha,'parent',r.parent); hold(r.parent,'on');
        
        switch(r.outline)
            case 1
                ph = plot(pnts,dens,'color',r.color, ...
                    'linew',r.linewidth,'linestyle',r.linestyle,'parent',r.parent); 
            otherwise
                ph = [];
        end

    case 'vertical'
        fh = fill(fr,pp,r.color, ...
            'EdgeColor',r.edgecolor,'FaceAlpha',r.facealpha,'parent',r.parent); hold(r.parent,'on');
        switch(r.outline)
            case 1
                ph = plot(dens,pnts,'color',r.color, ...
                    'linew',r.linewidth,'linestyle',r.linestyle,'parent',r.parent); 
            otherwise
                ph = [];
        end        

end



end