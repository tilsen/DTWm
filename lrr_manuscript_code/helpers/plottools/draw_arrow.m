function [arh] = draw_arrow(lh,loc,sz,edgecolor,facecolor)

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
PP = [0 2; -1 0; 0 0.25; 1 0];

%resize
PP = PP*sz;

%rotate
PP = PP*[cos(th) -sin(th); sin(th) cos(th)];

%displace
PP = PP+arxy;

arh = fill(PP(:,1),PP(:,2),facecolor,'EdgeColor',edgecolor);


end