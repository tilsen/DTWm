function [h] = plot_distance(r,q,map,varargin)

p = inputParser;

defaultgrid = false;
defaultgridcolor = [.5 .5 .5];
defaultmapcolor = [1 1 1];
defaultlinewidth = 2;
defaultmarker = 'none';
defaultlinestyle = '-';

q = q(:)';
r = r(:)';
if size(map,1)>size(map,2)
    map = map';
end   

addRequired(p,'r');
addRequired(p,'q');
addRequired(p,'map');
addParameter(p,'grid',defaultgrid);
addParameter(p,'gridcolor',defaultgridcolor);
addParameter(p,'mapcolor',defaultmapcolor);
addParameter(p,'linewidth',defaultlinewidth);
addParameter(p,'linestyle',defaultlinestyle);
addParameter(p,'marker',defaultmarker);

parse(p,r,q,map,varargin{:});

D = (q'-r).^2;

h.imh = imagesc(D); hold on;
axis tight;
if p.Results.grid
    h.gridlines = matrix_gridlines(gca,p.Results.gridcolor);
end

h.mh = plot(map(1,:),map(2,:),'linestyle',p.Results.linestyle, ...
    'marker',p.Results.marker,'linew',p.Results.linewidth,'color',p.Results.mapcolor);

set(gca,'YDir','normal');

end

