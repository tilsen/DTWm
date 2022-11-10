function [h] = plot_time_map(map,varargin)

p = inputParser;

def_linewidth = 2;
def_linetype = 'line';
def_colors = lines(2);
def_parent = nan;
def_markers = {'o','o'};
def_maplinecolor = [.5 .5 .5];
def_label_indices = true;
def_label_mappings = false;
def_fontsize = get(gca,'Fontsize');
def_yoffsets = [1 0];
def_yoffset_text = 0.25;

mapvalid = @(x)any(ismember(size(map),2)) & ismatrix(map);

addRequired(p,'map',mapvalid);
addParameter(p,'linewidth',def_linewidth);
addParameter(p,'linetype',def_linetype);
addParameter(p,'colors',def_colors);
addParameter(p,'markers',def_markers);
addParameter(p,'parent',def_parent);
addParameter(p,'maplinecolor',def_maplinecolor);
addParameter(p,'label_indices',def_label_indices);
addParameter(p,'label_mappings',def_label_mappings);
addParameter(p,'fontsize',def_fontsize);
addParameter(p,'yoffsets',def_yoffsets);
addParameter(p,'yoffset_text',def_yoffset_text);

parse(p,map,varargin{:});

r = p.Results;

%transpose if column vectors
if size(map,2)==2
    map = map';
end

s1 = map(1,:);
s2 = map(2,:);

s1u = unique(s1);
s2u = unique(s2);

if isnan(r.parent)
    r.parent = gca;
end

set(r.parent,'visible','off');
hold(r.parent,'on');

%%

th = [];
for i=1:size(map,2)
    
    %xx = [s1u(map(1,i)) s2u(map(2,i))];
    xx = [map(1,i) map(2,i)];
    yy = r.yoffsets;
    
    PP = [xx' yy'];
    
    switch(r.linetype)
        case 'arrow'
            ch(i) = arrow(PP(1,:),PP(2,:),...
                'length',8,'linewidth',r.linewidth,'tipangle',30,...
                'color',r.maplinecolor,'parent',r.parent); %#ok<*AGROW> 
            
        otherwise
            ch(i) = line(xx,yy,...
                'color',r.maplinecolor,'linew',r.linewidth,'parent',r.parent);
    end

    if r.label_mappings
        th(i) = text(mean(xx),mean(yy),num2str(i),'parent',r.parent,...
            'verti','mid','hori','right',...
            'fontsize',r.fontsize);
    end

end


yo = r.yoffset_text * -sign(diff(r.yoffsets));

th1= []; th2=[];

for i=1:length(s1u)
    lh1(i) = plot(s1u(i),r.yoffsets(1),'marker',r.markers{1},'MarkerFaceColor',r.colors(1,:),...
        'Color','k','parent',r.parent,'linew',1); 

    if r.label_indices
        th1(i) = text(s1u(i),r.yoffsets(1)+yo,num2str(s1u(i)),...
            'verti','top','hori','center','fontsize',r.fontsize);
    end
end

yo = -yo;
for i=1:length(s2u)
    lh2(i) = plot(s2u(i),r.yoffsets(2),'marker',r.markers{2},'MarkerFaceColor',r.colors(2,:),...
        'Color','k','parent',r.parent,'linew',1); 
    if r.label_indices
        th2(i) = text(s2u(i),r.yoffsets(2)+yo,num2str(s2u(i)),'verti','bot','hori','center');
    end    
end


h.ch = ch;
h.th = th;
h.lh1 = lh1;
h.lh2 = lh2;
h.th1 = th1;
h.th2 = th2;

axis(r.parent,'tight');
ylim(r.parent,ylim(r.parent)+[-1 1]);

set(r.parent,'XTick',[],'YTick',[]);


end