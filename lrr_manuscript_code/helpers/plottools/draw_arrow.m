function [arh] = draw_arrow(lh,varargin)

p = inputParser;

def_loc = 1;
def_sz = nan;
def_edgecolor = [0 0 0];
def_facecolor = [0 0 0];

addRequired(p,'lh',@(x)ishandle(x));
addOptional(p,'loc',def_loc,@(x)isscalar(x));
addOptional(p,'sz',def_sz,@(x)isscalar(x));
addParameter(p,'edgecolor',def_edgecolor);
addParameter(p,'facecolor',def_facecolor);

parse(p,lh,varargin{:});

r = p.Results;

loc = r.loc;
lh = r.lh;

if isnan(r.sz)
    sz = 0.01*diff(xlim(lh.Parent));
else
    sz = r.sz;
end

XY = [lh.XData' lh.YData'];

if numel(XY)>4 %circle
    arh = fill(nan,nan,[0 0 0]);
    return;
end

%vector
v = diff(XY);

%angle of vector
th = atan2(v(1),v(2));

%rescale vector
vs = v*loc;

%add to start point
arxy = XY(1,:) + vs;

%schematic arrow points (ToDo: implement more styles)
PP = [0 2; -1 0; 0 0.5; 1 0];

%resize
PP = PP*sz;

%rotate
PP = PP*[cos(th) -sin(th); sin(th) cos(th)];

%displace
PP = PP+arxy;

arh = fill(PP(:,1),PP(:,2),r.facecolor,'EdgeColor',r.edgecolor);


end