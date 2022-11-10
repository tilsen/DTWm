function [] = dtwm_validate()
%compares warping outputs to dtw-python for all step patterns

dbstop if error;

%load previously generated signals
fid = fopen(['.' filesep 'validation' filesep 'validation_signals'],'r');
for i=1:4, d{i} = str2num(fgetl(fid)); end; fclose(fid); %#ok<AGROW,ST2NM>

%%


%calculate warping for all step patterns
s = select_step_pattern; %no input -> return table of all patterns

%remove "rigid" b/c it will fail with different-length sequences
%unless open_end or open_begin is specified

s = s(~ismember(s.pattern,'rigid'),:);

%python-dtw does not recognize these step patterns (not sure why...)
s = s(~ismember(s.pattern,{'symmetricVelichkoZagoruyko' 'asymmetricItakura' 'asymmetricP0b'}),:);

tic
for i=1:height(s)
    s.warp{i} = dtwm(d{1},d{2},'step_pattern',s.pattern{i});
end
toc

fid = fopen(['.' filesep 'validation' filesep 'validation_dtwpython_output'],'r');
res = {};
while ~feof(fid)
    res{end+1} = fgets(fid);
end
fclose(fid);

for i=1:2:length(res)
   strs = regexp(res{i},'(\w+)(.*)','tokens'); strs = strs{1};
   R(i).pattern = strs{1};
   R(i).index1s = str2num(strs{2});
   strs = regexp(res{i+1},'(\w+)(.*)','tokens'); strs = strs{1};
   R(i).index2s = str2num(strs{2});
end
R = struct2table(R(1:2:end));

%compare python maps to dtwm maps
s.index1s = cellfun(@(c){c(1,:)},s.warp);
s.index2s = cellfun(@(c){c(2,:)},s.warp);

%adjust for python 0-indexing
R.index1s = cellfun(@(c){c+1},R.index1s);
R.index2s = cellfun(@(c){c+1},R.index2s);

for i=1:height(R)
    ix_s = find(ismember(s.pattern,R.pattern{i}));
    R.passed(i) = all(s.index1s{ix_s}==R.index1s{i}) && all(s.index2s{ix_s}==R.index2s{i});
end

disp(R(:,{'pattern' 'passed'}))

%% for validation of open end/begin
%to validate compare against values obtained in dtw_validate.py

patterns = {'symmetric2','asymmetric'};
opens = [false true];
for i=1:2
    for j=1:2
        [W{i,j}.map,W{i,j}.dist,W{i,j}.map_info] = ...
        dtwm(d{3},d{4},'open_begin',opens(i),'open_end',opens(j),'step_pattern',patterns{i});
    end
end


end





