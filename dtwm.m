function [wcurve,dist,gcm] = dtwm(x,varargin)
%This function "dtw-more" or "dtw-matlab" (dtwm) is a Matlab port written 
%by Sam Tilsen and Mark Tiede of Toni Giorgino's dtw package for python/R.
% See Giorgino (2009) Computing and Visualizing  Dynamic Time Warping 
% Alignments in R: The dtw Package.
%
%dtwm does more than built-in Matlab dtw: it allows for local slope 
% constraints, global windows, and partial matching
%
%Use: 
%     dtwm(x1,x2) where x1 and x2 are input signals
%        
% or  dtwm(D) where D is a distance matrix
% 
% Optional input parameters:
%
%   'dist_method':    'euclidean' (default); other methods not yet implemented
%
%   'step_pattern':   'symmetric2' (default); see documentation for alternatives
%
%   'window_type':    'none' (default), 'SakoeChiba', 'Itakura', 'SlantedBand'
%
%   'window_args':    {}; additional arguments for some windows. see
%                   documentation for examples.
%
%   'keep_internals': true (default) or false; keeps some information that might not
%                   be needed
%
%   'open_end':       true or false (default); allow for partial matching at end
% 
%   'open_begin':     true or false (default); allow for partial matching at begin
%
%   'use_mex':        true (default) or false; use compiled mex file for
%                   determining optimal warping path. 'false' will be very slow
%
% Outputs:
%
%   wcurve:     the warping curve
%
%   dist:       the distance after alignment
%
%   gcm:        additional information, including the global cost matrix



%{ 
The following is from the Python dtw.py function:

If you use this software in academic work, please cite:
#  * T. Giorgino. Computing and Visualizing Dynamic Time Warping
#    Alignments in R: The dtw Package. Journal of Statistical
#    Software, v. 31, Issue 7, p. 1 - 24, aug. 2009. ISSN
#    1548-7660. doi:10.18637/jss.v031.i07. http://www.jstatsoft.org/v31/i07/

Cost matrices (both input and output) have query elements arranged
row-wise (first index), and reference elements column-wise (second
index). They print according to the usual convention, with indexes
increasing down- and rightwards. Many DTW papers and tutorials show
matrices according to plot-like conventions, i_e. reference index
growing upwards. This may be confusing.

References
----------

1. Toni Giorgino. *Computing and Visualizing Dynamic Time Warping
   Alignments in R: The dtw Package.* Journal of Statistical Software,
   31(7), 1-24. http://www.jstatsoft.org/v31/i07/
2. Tormene, P.; Giorgino, T.; Quaglini, S. & Stefanelli, M. *Matching
   incomplete time series with dynamic time warping: an algorithm and an
   application to post-stroke rehabilitation.* Artif Intell Med, 2009,
   45, 11-34. http://dx.doi.org/10.1016/j.artmed.2008.11.007
3. Sakoe, H.; Chiba, S., *Dynamic programming algorithm optimization for
   spoken word recognition,* Acoustics, Speech, and Signal Processing,
   IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978.
   http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
4. Mori, A.; Uchida, S.; Kurazume, R.; Taniguchi, R.; Hasegawa, T. &
   Sakoe, H. *Early Recognition and Prediction of Gestures* Proc. 18th
   International Conference on Pattern Recognition ICPR 2006, 2006, 3,
   560-563
5. Sakoe, H. *Two-level DP-matching–A dynamic programming-based pattern
   matching algorithm for connected word recognition* Acoustics, Speech,
   and Signal Processing, IEEE Transactions on, 1979, 27, 588-595
6. Rabiner L, Rosenberg A, Levinson S (1978). *Considerations in dynamic
   time warping algorithms for discrete word recognition.* IEEE Trans.
   Acoust., Speech, Signal Process., 26(6), 575-582. ISSN 0096-3518.
7. Muller M. *Dynamic Time Warping* in *Information Retrieval for Music
   and Motion*. Springer Berlin Heidelberg; 2007. p. 69-84.
   http://link.springer.com/chapter/10.1007/978-3-540-74048-3_4

%}

p = inputParser;

default_dist_method = 'euclidean';
default_step_pattern = 'symmetric2';
default_window_type = '';
default_window_args = {};
default_keep_internals = true;
default_open_end = false;
default_open_begin = false;
default_use_mex = true;

validinput = @(c)(isnumeric(c) & numel(c)>1);
validinputy = @(c)(isnumeric(c) & numel(c)>1) | isempty(c);
addRequired(p,'x',validinput);
addOptional(p,'y',[],validinputy);
addParameter(p,'dist_method',default_dist_method);
addParameter(p,'step_pattern',default_step_pattern);
addParameter(p,'window_type',default_window_type);
addParameter(p,'window_args',default_window_args);
addParameter(p,'keep_internals',default_keep_internals);
addParameter(p,'open_end',default_open_end);
addParameter(p,'open_begin',default_open_begin);
addParameter(p,'use_mex',default_use_mex);

parse(p,x,varargin{:});
y = p.Results.y;

%get distance matrix (LM)
if isempty(y)
    if ~ismatrix(x)
        fprintf('error: first argument should be a matrix if no second argument is supplied\n'); return;
    end
    lm = x'; %transpose b/c treats rows as signal 1
else
    switch(p.Results.dist_method)
        case 'euclidean'
            if ndims(x)==1
                lm = sqrt((x'-y).^2);
            elseif ismatrix(x) && size(x,2)>size(x,1)
                lm = arrayfun(@(c){((x(c,:)'-y(c,:))).^2},(1:size(x,1)));
                lm = sqrt(sum(cat(3,lm{:}),3));
            else
                fprintf('error: signals must be vectors, or matrices with samples as second dimension\n');
                return;
            end
        otherwise
            fprintf('error: unknown distance method\n'); return;
    end
end

%get step pattern
if isstruct(p.Results.step_pattern)
    step_pattern = p.Results.step_pattern;
else
    step_pattern = select_step_pattern(p.Results.step_pattern);
end

%window function:
switch(lower(p.Results.window_type))
    case 'sakoechiba'
        if isempty(p.Results.window_args) || ~isscalar(p.Results.window_args)
            fprintf(['Error: Sakoe-Chiba requires a window size, e.g.\n' ...
                '\tdtwm(x1,x2,''window_type'',''sakoechiba'',''window_args'',20)\n']);
            return;
        end
        wfun=@(iw,jw,query_size,reference_size)sakoeChibaWindow(iw,jw,query_size,reference_size,p.Results.window_args);
    
    case 'itakura'

        wfun=@(iw,jw,query_size,reference_size)itakuraWindow(iw,jw,query_size,reference_size);
    
    case 'slantedband'

        wfun=@(iw,jw,query_size,reference_size)slantedBandWindow(iw,jw,query_size,reference_size,p.Results.window_args);
    
    otherwise
        if ~isempty(p.Results.window_type), fprintf('unrecognized window function\n'); return; end
        wfun=[];
end

norm = step_pattern.norm;
[n,m] = size(lm);

%prep inputs
if p.Results.open_begin
    if ~strcmp(norm,'N')
        fprintf('error: open begin requires step patterns with "N" normalization\n'); return;
    end
    lm = [zeros(1,size(lm,2)); lm]; %prepend null row
    np = n+1;
    precm = nan(size(lm));
    precm(1,:) = 0;
else
    precm = [];
    np = n;
end

%calculate cost matrix
switch(p.Results.use_mex)
    case true
        wm = true(size(lm));
        [n,m] = size(lm);
        if ~isempty(wfun)
            [ii,jj] = meshgrid(1:n,1:m);
            wm = wfun(ii,jj,n,m)';
        end
        
        pn = int32(step_pattern.mx(:,1)')-1;
        di = int32(step_pattern.mx(:,2)');
        dj = int32(step_pattern.mx(:,3)');
        sc = double(step_pattern.mx(:,4)');
        
        precm = nan(size(lm));
        precm(1,1) = 0;
        
        clist = double(nan(numel(unique(pn)),1))';
        sm = int32(nan(size(precm)));
        
        [gcm.costMatrix,gcm.directionMatrix] = dtwm_costmatrix(precm,pn,di,dj,sc,lm,clist,sm,wm);
        gcm.stepPattern = step_pattern;
        
    otherwise
        gcm = globalCostMatrix(lm,step_pattern,wfun,precm);
end

gcm.N = n;
gcm.M = m;
gcm.openEnd = p.Results.open_end;
gcm.openBegin = p.Results.open_begin;
gcm.windowFunction = wfun;
gcm.windowArgs = p.Results.window_args;

lastrow = gcm.costMatrix(end,:);

switch(norm)
    case 'N+M'
        lastrow = lastrow ./ (n + (1:m));
    case 'N'
        lastrow = lastrow / n;
    case 'M'
        lastrow = lastrow ./ (1:m);
end

gcm.jmin = m;

if p.Results.open_end
    if isempty(norm)
        fprintf('error: open end requires requires normalizable step patterns\n'); return;
    end
    [~,gcm.jmin] = min(lastrow);
end

gcm.distance = gcm.costMatrix(end,gcm.jmin);

if isnan(gcm.distance)
    fprintf('error: no warping path found compatible with the local constraints\n'); 
    wcurve = [];
    dist = nan;
    return;
end

if ~isempty(step_pattern.norm)
    gcm.normalizedDistance = lastrow(gcm.jmin);
else
    gcm.normalizedDistance = nan;
end

mapping = backtrack(gcm);

gcm.index1 = mapping.index1;
gcm.index2 = mapping.index2;
gcm.index1s = mapping.index1s;
gcm.index2s = mapping.index2s;

if p.Results.open_begin
    gcm.index1 = gcm.index1(2:end) - 1; 
    gcm.index1s = gcm.index1s(2:end) - 1;
    gcm.index2 = gcm.index2(2:end);
    gcm.index2s = gcm.index2s(2:end);
    lm = lm(2:end,:);
    gcm.costMatrix = gcm.costMatrix(2:end,:);
    gcm.directionMatrix = gcm.directionMatrix(2:end,:);
end

if ~p.Results.keep_internals
    gcm = rmfield(gcm,{'costMatrix','directionMatrix'});
else
    gcm.localCostMatrix = lm;
end

%outputs

wcurve = [gcm.index1s; gcm.index2s];
dist = gcm.normalizedDistance;

% adjust mapping for open end/begin:
if gcm.openEnd && ~gcm.openBegin
    %jmin holds the size of the prefixed matched, the index j of the
    %query
    wcurve = wcurve(:,1:find(wcurve(2,:)==gcm.jmin,1,'first'));

elseif ~gcm.openEnd && gcm.openBegin
    wcurve = wcurve(:,find(wcurve(1,:)>0,1,'first'):end);
    wcurve = wcurve(:,find(wcurve(2,:)==wcurve(2,1),1,'last'):end);

elseif gcm.openEnd && gcm.openBegin

    wcurve = wcurve(:,1:find(wcurve(2,:)==gcm.jmin,1,'first'));
    wcurve = wcurve(:,find(wcurve(1,:)>0,1,'first'):end);
    wcurve = wcurve(:,find(wcurve(2,:)==wcurve(2,1),1,'last'):end);
end

end

%% backtrack
function [mapping] = backtrack(gcm)

i = gcm.N;
j = gcm.M;

dir = gcm.stepPattern.mx;

%drop null deltas
dir = dir(dir(:,2)~=0 | dir(:,3)~=0,:);

npat = numel(unique(dir(:,1)));

stepsCache = {};

for q=1:npat
    tmp = dir(dir(:,1)==q,:);
    stepsCache{q} = flipud(tmp(:,[2 3])); %#ok<AGROW> 
end

%mapping lists:
iis = i;
iix = i;
jjs = j;
jjx = j;
ss = [];

while 1
    if i==1 && j==1, break; end
    
    %direction taken, 1-based
    s = gcm.directionMatrix(i,j);
    
    if isnan(s), break; end
    
    ss = [s ss];
    
    steps = stepsCache{s};
    
    for k=1:size(steps,1)
        iix = [i - steps(k,1) iix];
        jjx = [j - steps(k,2) jjx];
    end
    
    i = i - steps(k,1);
    j = j - steps(k,2);
    
    iis = [i iis];
    jjs = [j jjs];
    
end

mapping.index1 = iix;
mapping.index2 = jjx;
mapping.index1s = iis;
mapping.index2s = jjs;

end

%% cost matrix
function [gcm] = globalCostMatrix(lm,step_pattern,window_function,seed)

wm = true(size(lm));
[n,m] = size(wm);

if ~isempty(window_function)
    if ~isempty(wfun)
        [ii,jj] = meshgrid(1:n,1:m);
        wm = window_function(ii,jj,n,m)';
    end
end

nsteps = size(step_pattern.mx,1);

pn = step_pattern.mx(:,1); %pattern id
di = step_pattern.mx(:,2); %delta i
dj = step_pattern.mx(:,3); %delta j
sc = step_pattern.mx(:,4); %step cost

%cost matrix
if ~isempty(seed)
    cm = seed;
else
    cm = nan(size(lm));
    cm(1,1) = lm(1,1);
end

%steps/direction matrix
dm = nan(size(cm));

[m,n] = size(lm);

npats = numel(unique(pn));

for j=1:n %loop over columns
    for i=1:m %loop over rows
        costs = nan(npats,1);
        for s=1:nsteps
            p = pn(s);
            ii = i-di(s);
            jj = j-dj(s);
            if (ii>=1 && jj>=1)
                cc = sc(s);
                if cc==-1
                    costs(p) = cm(ii,jj);
                else
                    costs(p) = costs(p) + cc*lm(ii,jj); %relies on NaN propagation here
                end
            end
        end
        
        if ~all(isnan(costs))
            [~,minc] = min(costs);
            cm(i,j) = costs(minc);
            dm(i,j) = minc;
        end
    end
end

gcm.costMatrix = cm;
gcm.directionMatrix = dm;
gcm.stepPattern = step_pattern;

end

%% window functions:
function [out] = sakoeChibaWindow(iw,jw,~,~,window_size)
out = abs(jw-iw)<=window_size;
end

function [out] = itakuraWindow(iw,jw,query_size,reference_size,~)
n = query_size;
m = reference_size;
jw = jw-1;
iw = iw-1;
out = (jw < 2*iw) & (iw<= 2*jw) & (iw>=n -1 - 2*(m-jw)) & (jw > m - 1 - 2*(n-iw));
end

function [out] = slantedBandWindow(iw,jw,query_size,reference_size,window_size)
diagj = (iw * reference_size / query_size);
out = abs(jw-diagj) <= window_size;
end

