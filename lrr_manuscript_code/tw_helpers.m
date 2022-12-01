function [h] = tw_helpers()


cd(fileparts(mfilename('fullpath')));
addpath(genpath('.'));

%define data directory (should contain onsetcoda_opendata.mat)
switch(ispc)
    case 1
        h.data_dir = 'M:\Data\timewarping_opendata\';        
    otherwise
        h.data_dir = '/home/tilsen/Data/timewarping_opendata/';      
        addpath('/home/tilsen/Projects/toolboxes/st_tools/dtw_more/'); 
end

if ~exist(h.data_dir,'dir')
    fprintf('Data directory: %s does not exist\n',h.data_dir); return;
end

%define optional toolbox paths here:
toolboxes = {
    'M:\Projects\toolboxes\arrow\' 'https://www.mathworks.com/matlabcentral/fileexchange/278-arrow'
    'M:\Projects\toolboxes\drawbrace\' 'https://www.mathworks.com/matlabcentral/fileexchange/38716-curly-brace-annotation'
    'M:\Projects\toolboxes\st_tools\dtw_more\' 'https://github.com/tilsen/DTWm.git'
    };

for i=1:size(toolboxes,1)
    if ~exist(toolboxes{i,1},'dir')
        %fprintf('\nWARNING: The following toolbox may be required:\n%s\n',toolboxes{i,2});
    else
        addpath(toolboxes{i,1});
    end
end

h.datasets_dir = [h.data_dir 'datasets' filesep];
h.lrr_dir = [h.data_dir 'LRR' filesep];
h.gam_dir = [h.data_dir 'GAM' filesep];

%%
h.fs = [36 24 20 16]; %fontsizes
h.ticks = {'tickdir','out','ticklen',[0.0030 0.0030]}; %default axes ticks

%%

h.mathsym = @(str)['{\fontname{CMU Serif}{\it' str '}}'];

h.printfig = @(str)print('-dpng','-r350',['.' filesep 'images' filesep str '.png']);

end
