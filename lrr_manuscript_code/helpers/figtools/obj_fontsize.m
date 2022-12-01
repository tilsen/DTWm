function [] = obj_fontsize(obj,varargin)

object_types = varargin(1:2:end);
font_sizes = cell2mat(varargin(2:2:end));

ch = allchild(obj);
%ch = vertcat(ch{:});
ax = ch(ismember(get(ch,'type'),'axes'));

CH = allchild(ax);
ch = [ch; vertcat(CH{:})];

for i=1:length(object_types)

    switch(object_types{i})
        case 'label'
            for j=1:length(ax)
                ax(j).XLabel.FontSize = font_sizes(i);
                ax(j).YLabel.FontSize = font_sizes(i);
                ax(j).ZLabel.FontSize = font_sizes(i);
            end

        case 'xlabel'
            for j=1:length(ax)
                ax(j).XLabel.FontSize = font_sizes(i);
            end            

        case 'ylabel'
            for j=1:length(ax)
                ax(j).YLabel.FontSize = font_sizes(i);
            end            

        otherwise
            set(ch(ismember(get(ch,'type'),object_types{i})),'Fontsize',font_sizes(i));
    end

end

end