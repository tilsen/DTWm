function [slopes,singular_warn] = lrr(map,varargin)
%lrr computes the local rate of the first component of the map relative to
%the second component.
%
%Input:
%
%   first row/col of map is usually the target/reference/abscissa/domain of dtw
%   warping curve
%
%   second row/col of map is usually comparison/test/ordinate/codomain of dtw
%   warping curve
%
%Output:
%
%   slopes: slope of change in test warping curve index to 
%                 change in ref warping curve index
%
% equivalent to:
%
%   local relative rate (lrr) of ref/target index relative to
%   text/comparison index

p = inputParser;

def_method = 'regression';
def_window = 'triangular';   % triangular, gaussian;
def_winsize = 0.05;          % proportion of length if <1
def_warn_opt = 'catch';      % disable, catch (records if local regression matrix is nearly singular) 

mapvalid = @(x)any(ismember(size(map),2)) & ismatrix(map);

addRequired(p,'map',mapvalid);
addParameter(p,'method',def_method);
addParameter(p,'window',def_window);
addParameter(p,'winsize',def_winsize);
addParameter(p,'warn_opt',def_warn_opt);

parse(p,map,varargin{:});

r = p.Results;

%transpose if column vectors
if size(map,2)==2
    map = map';
end

s1 = map(1,:);
%s2 = map(2,:);

% time points (indices)
uci1 = min(s1):max(s1);    %unique consecutive indices vector
L = length(uci1);          %length of unique consecutive indices vector
N = length(s1);            %length of map

%determine winsize if not supplied
if r.winsize<1
    r.winsize = ceil(r.winsize*L);
end

%warnings
warnstrs = {'MATLAB:nearlySingularMatrix','MATLAB:singularMatrix'};
switch(r.warn_opt)
    case 'disable'
        warn_opt = 'off';
    case 'catch'
        warn_opt = 'error';
    otherwise
        warn_opt = 'on';
end
for i=1:length(warnstrs)
    warning(warn_opt,warnstrs{i});
end

%% calculate weights
switch(r.window)
    case 'triangular'
        wfcn = @(t,tc)max(0,1-abs((t-tc)/(r.winsize/2)));
    case 'gaussian'
        wfcn = @(t,tc)exp(-((t-tc).^2)/(2*(r.winsize/2)^2));
    case 'rectangular'
        wfcn = @(t,tc)abs((t-tc))<=(r.winsize/2);
end

wx = speye(N);
w = cell(1,L);
for i=1:L
    w{i} = wx.*wfcn(s1,uci1(i));
end

%% calculate slopes
coeffs = nan(2,L);
singular_warn = false(1,L);

m = map-mean(map,2); %re-center (subtract means)
m = m/max(abs(m(:))); %re-scale (divide by maximal abs value)

X = [ones(1,N); m(1,:)]';
Y = m(2,:)';

switch(r.method)
    case 'regression'
        for i=1:L
            try
                coeffs(:,i) = (X'*w{i}*X)\(X'*w{i}*Y);

            catch
                singular_warn(i) = true;

                for j=1:length(warnstrs)
                    warning('off',warnstrs{j});
                end
                coeffs(:,i) = (X'*w{i}*X)\(X'*w{i}*Y);
                for j=1:length(warnstrs)
                    warning('error',warnstrs{j});
                end
            end
        end

end

slopes = coeffs(2,:);

end


