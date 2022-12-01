function [obj] = graph_plot(obj,varargin)
%does what GraphPlot does but puts individual node, text, and line objects
%in tables

p = inputParser;

def_NodeSize = 10;
def_NodeFaceColor = [1 1 1];
def_NodeLineColor = [0 0 0];
def_ArrowPosition = 0.5;
def_ArrowSize = 0.02; %data units
def_EdgeWidth = 1;
def_EdgeColor = [0 0 0];
def_SelfLoopRadius = 0.1;
def_NodeLabelFontSize = 12;
def_NodeLabelPosition = 0.5;
def_EdgeStyle = '-';
def_ArrowEdgeColor = [0 0 0];
def_ArrowFaceColor = [0 0 0];

addRequired(p,'obj');
addParameter(p,'XData',[]);
addParameter(p,'YData',[]);
addParameter(p,'NodeSize',def_NodeSize);
addParameter(p,'NodeFaceColor',def_NodeFaceColor);
addParameter(p,'NodeLineColor',def_NodeLineColor);
addParameter(p,'NodeLabelPosition',def_NodeLabelPosition);
addParameter(p,'NodeLabelFontSize',def_NodeLabelFontSize);
%addParameter(p,'NodeLabelHorizontalAlignment',def_NodeLabelHorizontalAlignment);

addParameter(p,'EdgeWidth',def_EdgeWidth);
addParameter(p,'EdgeStyle',def_EdgeStyle);
addParameter(p,'EdgeColor',def_EdgeColor);
addParameter(p,'SelfLoopRadius',def_SelfLoopRadius);

addParameter(p,'ArrowPosition',def_ArrowPosition);
addParameter(p,'ArrowSize',def_ArrowSize);
addParameter(p,'ArrowEdgeColor',def_ArrowEdgeColor);
addParameter(p,'ArrowFaceColor',def_ArrowFaceColor);

parse(p,obj,varargin{:});

r = p.Results;

% if Xdata and Ydata not provided, use GraphPlot layout
if isempty(r.XData) || isempty(r.YData)
    h = plot(obj);
    r.XData = h.XData;
    r.YData = h.YData;
    delete(h);
    clear(h);
else
    if numel(r.XData)~=height(obj.Nodes) || numel(r.XData)~=numel(r.YData)
        fprintf('error: XData and YData must be specified for each node\n');
        return;
    end
end



getxy = @(e)[r.XData(ismember(obj.Nodes.Name,e)) r.YData(ismember(obj.Nodes.Name,e))];
linexy = @(e)[getxy(e(1)); getxy(e(2))];

%draw edges
for i=1:height(obj.Edges)

    xy = linexy(obj.Edges.EndNodes(i,:));

    if all(diff(xy)==0)
        th = linspace(0,2*pi,100);
        cxy = r.SelfLoopRadius*[cos(th)' sin(th)'];
        xy = mean(xy)+cxy+r.SelfLoopRadius*[-1 1];
    end

    obj.Edges.lh(i) = plot(xy(:,1),xy(:,2), ...
        'Color',r.EdgeColor, ...
        'Linewidth',r.EdgeWidth, ...
        'LineStyle',r.EdgeStyle); hold on;

    obj.Edges.arh(i) = draw_arrow(obj.Edges.lh(i), ...
        r.ArrowPosition, ...
        r.ArrowSize, ...
        r.ArrowEdgeColor,...
        r.ArrowFaceColor);
end


%draw nodes
for i=1:height(obj.Nodes)
    obj.Nodes.nh(i) = plot(r.XData(i),r.YData(i),'o', ...
        'color',r.NodeLineColor, ...
        'MarkerFaceColor',r.NodeFaceColor,...
        'MarkerSize',r.NodeSize); hold on;

    obj.Nodes.th(i) = text(r.XData(i),r.YData(i),obj.Nodes.Name{i}, ...
        'hori','center','fontsize',r.NodeLabelFontSize);

    
end

% add edgelabels
axis tight;
axis equal;

for i=1:height(obj.Edges)
    xy = [obj.Edges.lh(i).XData' obj.Edges.lh(i).YData'];

    if numel(xy)>4
        xy = mean(xy);
        P = 0;

    else
        aa = diff(xy);

        P = rad2deg(atan2(aa(2),aa(1)));
        P(P>90) = P(P>90)-180;
        P(P<-90) = P(P<-90)+180;

        xy = xy(1,:)+r.NodeLabelPosition*diff(xy);
    end

    obj.Edges.th(i) = text(xy(1),xy(2),tf(obj.Edges.Weight(i)), ...
        'FontSize',12, ...
        'Verti','bot',...
        'hori','center',...
        'Rotation',P);

end

%collect handles
for i=1:height(obj.Edges)
    obj.Edges.handles(i,:) = [obj.Edges.lh(i) obj.Edges.arh(i) obj.Edges.th(i)];
end

end

%format text
function [str] = tf(x)
str = num2str(x);

end

