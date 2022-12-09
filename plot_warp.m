function [h] = plot_warp(q,r,map,varargin)

p = inputParser;

defaultdirection = 'bidirectional';
defaultcolormap = 'auto';
defaultyoffset = 'auto';
defaultmaparrow = 'none';
defaultmode = 'unwarped';

q = q(:)';
r = r(:)';
if size(map,1)>size(map,2)
    map = map';
end

addRequired(p,'q');
addRequired(p,'r');
addRequired(p,'map');
addParameter(p,'direction',defaultdirection);
addParameter(p,'colormap',defaultcolormap);
addParameter(p,'yoffset',defaultyoffset);
addParameter(p,'mapline',defaultmaparrow);
addParameter(p,'mode',defaultmode);

parse(p,q,r,map,varargin{:});

switch(p.Results.colormap)
    case 'auto'
        rng = minmax([q r]);
        cmap = jet(1000);
        qn = q-min(rng); qn = qn/diff(rng);
        rn = r-min(rng); rn = rn/diff(rng);
        qcolors = cmap(1+floor(qn*999),:);
        rcolors = cmap(1+floor(rn*999),:);
    otherwise
        
end


switch(p.Results.yoffset)
    case 'auto'
        yo = max(r);
    otherwise
        yo = p.Results.yoffset;
end

%check for arrow toolbox:
artoolbox = which('arrow');
if isempty(artoolbox)
    res.mapline = 'none';
else
    res = p.Results;
end

switch(p.Results.mode)
    case 'warped'
        q = q(map(1,:));
        r = r(map(2,:));
        
        h.q_line = plot(1:length(q),q,'-','color','k','linew',3); hold on;
        h.r_line = plot(1:length(r),r-yo,'-','color','k','linew',3); hold on;
        for i=1:length(q)
            h.q(i) = plot(i,q(i),'ko','markerfacecolor',qcolors(map(1,i),:),'markersize',10); hold on;
        end
        for i=1:length(r)
            h.r(i) = plot(i,r(i)-yo,'ko','markerfacecolor',rcolors(map(2,i),:),'markersize',10); hold on;
        end        
        
    otherwise
        
        h.q_line = plot(1:length(q),q,'-','color','k','linew',3); hold on;
        h.r_line = plot(1:length(r),r-yo,'-','color','k','linew',3); hold on;
        for i=1:length(q)
            h.q(i) = plot(i,q(i),'ko','markerfacecolor',qcolors(i,:),'markersize',10); hold on;
        end
        for i=1:length(r)
            h.r(i) = plot(i,r(i)-yo,'ko','markerfacecolor',rcolors(i,:),'markersize',10); hold on;
        end
        
end

axis tight;
xlim(xlim);
ylim(ylim);

for i=1:size(map,2)

    switch(p.Results.mode)
        case 'warped'
            PP = [  h.q(i).XData  h.q(i).YData;...
                    h.r(i).XData  h.r(i).YData];            
        otherwise
            PP = [  h.q(map(1,i)).XData  h.q(map(1,i)).YData;...
                h.r(map(2,i)).XData  h.r(map(2,i)).YData];
    end
    
    PP = scale_connection(PP,0.90);
    
    switch(res.mapline)
        case 'none'
            h.mh(i) = line(PP(:,1),PP(:,2),'linewidth',2,'color','k');
            
        case 'bidirectional'
            h.mh(i,1) = arrow(PP(1,:),PP(2,:),'length',8,'linewidth',2,'tipangle',30);
            h.mh(i,2) = arrow(PP(2,:),PP(1,:),'length',8,'linewidth',2,'tipangle',30);
        
        case 'unidirectional'
            h.mh(i,1) = arrow(PP(1,:),PP(2,:),'length',8,'linewidth',2,'tipangle',30);
                    
    end
    
end

end

function [Ca] = scale_connection(C,f)

%stretches or contracts a line
muC = mean(C);

%center line at origin
Cn = C-repmat(muC,2,1);

%scale
Cn = Cn*f;

%restore position
Ca = Cn+repmat(muC,2,1);


end

