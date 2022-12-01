function [ph,axm] = plot_lrr_array(X,colors)

lw = 2;

Ns = size(X,1);

axm = stfig_subaxpos(gca,[size(X)],[0 0 0 0 0 0]);

axm = reshape(axm,Ns,[])';

L = cellfun('length',X);
maxL = max(L(:));
maxX = max([X{:}]);
minX = min([X{:}]);

for r=1:size(X,1)
    for c=1:size(X,2)
        plot(1:maxL,ones(1,maxL),'-','color',[.5 .5 .5],'parent',axm(r,c)); hold(axm(r,c),'on');
        ph(r,c) = plot(X{r,c},'-','color',colors(c,:),'linew',lw,'parent',axm(r,c));
    end
end

set(axm,'YLim',[minX maxX],'XLim',[0 maxL],'YTick',[],'XTick',[]);
axrescale(axm,0.05,0.05);

end