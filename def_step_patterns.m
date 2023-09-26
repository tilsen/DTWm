function [s] = def_step_patterns()


% ##
% ## Various step patterns, defined as internal variables
% ##
% ## First column: enumerates step patterns.
% ## Second   	 step in query index
% ## Third	 step in reference index
% ## Fourth	 weight if positive, or -1 if starting point
% ##
% ## For \cite{} see dtw.bib in the package
% ##

%
% ## Widely-known variants
%
% ## White-Neely symmetric (default)
% ## aka Quasi-symmetric \cite{White1976}
% ## normalization: no (N+M?)
s(1).pattern = 'symmetric1';
s(end).mx = [
    1, 1, 1, -1;
    1, 0, 0, 1;
    2, 0, 1, -1;
    2, 0, 0, 1;
    3, 1, 0, -1;
    3, 0, 0, 1];

% ## Normal symmetric
% ## normalization: N+M
s(end+1).pattern = 'symmetric2'; 
s(end).mx = [...
    1, 1, 1, -1;
    1, 0, 0, 2;
    2, 0, 1, -1;
    2, 0, 0, 1;
    3, 1, 0, -1;
    3, 0, 0, 1];
s(end).norm = "N+M";

% ## classic asymmetric pattern: max slope 2, min slope 0
% ## normalization: N
s(end+1).pattern = 'asymmetric';
s(end).mx = [...
    1, 1, 0, -1;
    1, 0, 0, 1;
    2, 1, 1, -1;
    2, 0, 0, 1;
    3, 1, 2, -1;
    3, 0, 0, 1];
s(end).norm = "N";

% # % \item{\code{symmetricVelichkoZagoruyko}}{symmetric, reproduced from %
% # [Sakoe1978]. Use distance matrix \code{1-d}}
% #
%
% ## normalization: max[N,M]
% ## note: local distance matrix is 1-d
% ## \cite{Velichko}
s(end+1).pattern = 'symmetricVelichkoZagoruyko';
s(end).mx = [...
    1, 0, 1, -1;
    2, 1, 1, -1;
    2, 0, 0, -1.001;
    3, 1, 0, -1];

% # % \item{\code{asymmetricItakura}}{asymmetric, slope contrained 0.5 -- 2
% # from reference [Itakura1975]. This is the recursive definition % that
% # generates the Itakura parallelogram; }
% #
%
% ## Itakura slope-limited asymmetric \cite{Itakura1975}
% ## Max slope: 2; min slope: 1/2
% ## normalization: N
s(end+1).pattern = 'asymmetricItakura';
s(end).mx = [...
    1, 1, 2, -1;
    1, 0, 0, 1;
    2, 1, 1, -1;
    2, 0, 0, 1;
    3, 2, 1, -1;
    3, 1, 0, 1;
    3, 0, 0, 1;
    4, 2, 2, -1;
    4, 1, 0, 1;
    4, 0, 0, 1];

% #############################
% ## Slope-limited versions
% ##
% ## Taken from Table I, page 47 of "Dynamic programming algorithm
% ## optimization for spoken word recognition," Acoustics, Speech, and
% ## Signal Processing, vol.26, no.1, pp. 43-49, Feb 1978 URL:
% ## http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
% ##
% ## Mostly unchecked
%
%
% ## Row P=0
s(end+1).pattern = 'symmetricP0';
s(end).mx = [...
    1, 1, 1, -1;
    1, 0, 0, 2;
    2, 0, 1, -1;
    2, 0, 0, 1;
    3, 1, 0, -1;
    3, 0, 0, 1];
s(end).norm = "N+M";
%
% ## normalization: N ?
s(end+1).pattern = 'asymmetricP0';
s(end).mx = [...
    1, 0, 1, -1;
    1, 0, 0, 0;
    2, 1, 1, -1;
    2, 0, 0, 1;
    3, 1, 0, -1;
    3, 0, 0, 1];
s(end).norm = "N";

%## alternative implementation
s(end+1).pattern = 'asymmetricP0b';
s(end).mx = [...
    1, 0, 1, -1;
    2, 1, 1, -1;
    2, 0, 0, 1;
    3, 1, 0, -1;
    3, 0, 0, 1];
s(end).norm = "N";

%## Row P=1/2
s(end+1).pattern = 'symmetricP05';
s(end).mx = [...
    1, 1, 3, -1;
    1, 0, 2, 2;
    1, 0, 1, 1;
    1, 0, 0, 1;
    2, 1, 2, -1;
    2, 0, 1, 2;
    2, 0, 0, 1;
    3, 1, 1, -1;
    3, 0, 0, 2;
    4, 2, 1, -1;
    4, 1, 0, 2;
    4, 0, 0, 1;
    5, 3, 1, -1;
    5, 2, 0, 2;
    5, 1, 0, 1;
    5, 0, 0, 1];
s(end).norm = "N+M";

s(end+1).pattern = 'asymmetricP05';
s(end).mx = [...
    1, 1, 3, -1;
    1, 0, 2, 1 / 3;
    1, 0, 1, 1 / 3;
    1, 0, 0, 1 / 3;
    2, 1, 2, -1;
    2, 0, 1, .5;
    2, 0, 0, .5;
    3, 1, 1, -1;
    3, 0, 0, 1;
    4, 2, 1, -1;
    4, 1, 0, 1;
    4, 0, 0, 1;
    5, 3, 1, -1;
    5, 2, 0, 1;
    5, 1, 0, 1;
    5, 0, 0, 1];
s(end).norm = "N";

% ## Row P=1
% ## Implementation of Sakoe's P=1, Symmetric algorithm

s(end+1).pattern = 'symmetricP1';
s(end).mx = [...
    1, 1, 2, -1; ...  # First branch: g(i-1,j-2)+
    1, 0, 1, 2; ...  # + 2d(i  ,j-1)
    1, 0, 0, 1; ...  # +  d(i  ,j)
    2, 1, 1, -1; ...  # Second branch: g(i-1,j-1)+
    2, 0, 0, 2; ...  # +2d(i,  j)
    3, 2, 1, -1; ...  # Third branch: g(i-2,j-1)+
    3, 1, 0, 2; ...  # + 2d(i-1,j)
    3, 0, 0, 1]; %  # +  d(  i,j)
s(end).norm = "N+M";

s(end+1).pattern = 'asymmetricP1';
s(end).mx = [...
    1, 1, 2, -1; ...
    1, 0, 1, .5;
    1, 0, 0, .5;
    2, 1, 1, -1;
    2, 0, 0, 1;
    3, 2, 1, -1;
    3, 1, 0, 1;
    3, 0, 0, 1];
s(end).norm = "N";

% ## Row P=2
s(end+1).pattern = 'symmetricP2';
s(end).mx = reshape([
    1, 2, 3, -1,
    1, 1, 2, 2,
    1, 0, 1, 2,
    1, 0, 0, 1,
    2, 1, 1, -1,
    2, 0, 0, 2,
    3, 3, 2, -1,
    3, 2, 1, 2,
    3, 1, 0, 2,
    3, 0, 0, 1],[],4);
s(end).norm = "N+M";

s(end+1).pattern = 'asymmetricP2';
s(end).mx = reshape([
    1, 2, 3, -1,
    1, 1, 2, 2 / 3,
    1, 0, 1, 2 / 3,
    1, 0, 0, 2 / 3,
    2, 1, 1, -1,
    2, 0, 0, 1,
    3, 3, 2, -1,
    3, 2, 1, 1,
    3, 1, 0, 1,
    3, 0, 0, 1],[],4);
s(end).norm = "N";

% ################################
% ## Taken from Table III, page 49.
% ## Four varieties of DP-algorithm compared
%
% ## 1st row:  asymmetric
%
% ## 2nd row:  symmetricVelichkoZagoruyko
%
% ## 3rd row:  symmetric1
%
% ## 4th row:  asymmetricItakura
%
%
% #############################
% ## Classified according to Rabiner
% ##
% ## Taken from chapter 2, Myers' thesis [4]. Letter is
% ## the weighting function:
% ##
% ##      rule       norm   unbiased
% ##   a  min step   ~N     NO
% ##   b  max step   ~N     NO
% ##   c  x step     N      YES
% ##   d  x+y step   N+M    YES
% ##
% ## Mostly unchecked
%
% # R-Myers     R-Juang
% # type I      type II
% # type II     type III
% # type III    type IV
% # type IV     type VII


s(end+1).pattern = 'typeIa';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, 1,
    1, 0, 0, 0,
    2, 1, 1, -1,
    2, 0, 0, 1,
    3, 1, 2, -1,
    3, 0, 1, 1,
    3, 0, 0, 0],[],4);

s(end+1).pattern = 'typeIb';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, 1,
    1, 0, 0, 1,
    2, 1, 1, -1,
    2, 0, 0, 1,
    3, 1, 2, -1,
    3, 0, 1, 1,
    3, 0, 0, 1],[],4);

s(end+1).pattern = 'typeIc';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, 1,
    1, 0, 0, 1,
    2, 1, 1, -1,
    2, 0, 0, 1,
    3, 1, 2, -1,
    3, 0, 1, 1,
    3, 0, 0, 0],[],4);
s(end).norm = "N";

s(end+1).pattern = 'typeId';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, 2,
    1, 0, 0, 1,
    2, 1, 1, -1,
    2, 0, 0, 2,
    3, 1, 2, -1,
    3, 0, 1, 2,
    3, 0, 0, 1],[],4);
s(end).norm = "N+M";

% ## ----------
% ## smoothed variants of above

s(end+1).pattern = 'typeIas';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, .5,
    1, 0, 0, .5,
    2, 1, 1, -1,
    2, 0, 0, 1,
    3, 1, 2, -1,
    3, 0, 1, .5,
    3, 0, 0, .5],[],4);

s(end+1).pattern = 'typeIbs';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, 1,
    1, 0, 0, 1,
    2, 1, 1, -1,
    2, 0, 0, 1,
    3, 1, 2, -1,
    3, 0, 1, 1,
    3, 0, 0, 1],[],4);

s(end+1).pattern = 'typeIcs';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, 1,
    1, 0, 0, 1,
    2, 1, 1, -1,
    2, 0, 0, 1,
    3, 1, 2, -1,
    3, 0, 1, .5,
    3, 0, 0, .5],[],4);
s(end).norm = "N";

s(end+1).pattern = 'typeIds';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, 1.5,
    1, 0, 0, 1.5,
    2, 1, 1, -1,
    2, 0, 0, 2,
    3, 1, 2, -1,
    3, 0, 1, 1.5,
    3, 0, 0, 1.5],[],4);
s(end).norm =  "N+M";

%## ----------

s(end+1).pattern = 'typeIIa';
s(end).mx = reshape([
    1, 1, 1, -1,
    1, 0, 0, 1,
    2, 1, 2, -1,
    2, 0, 0, 1,
    3, 2, 1, -1,
    3, 0, 0, 1],[],4);

s(end+1).pattern = 'typeIIb';
s(end).mx = reshape([
    1, 1, 1, -1,
    1, 0, 0, 1,
    2, 1, 2, -1,
    2, 0, 0, 2,
    3, 2, 1, -1,
    3, 0, 0, 2],[],4);

s(end+1).pattern = 'typeIIc';
s(end).mx = reshape([
    1, 1, 1, -1,
    1, 0, 0, 1,
    2, 1, 2, -1,
    2, 0, 0, 1,
    3, 2, 1, -1,
    3, 0, 0, 2],[],4);
s(end).norm = "N";

s(end+1).pattern = 'typeIId';
s(end).mx = reshape([
    1, 1, 1, -1,
    1, 0, 0, 2,
    2, 1, 2, -1,
    2, 0, 0, 3,
    3, 2, 1, -1,
    3, 0, 0, 3],[],4);
s(end).norm = "N+M";

% ## ----------
%
% ## Rabiner [3] discusses why this is not equivalent to Itakura's

s(end+1).pattern = 'typeIIIc';
s(end).mx = reshape([
    1, 1, 2, -1,
    1, 0, 0, 1,
    2, 1, 1, -1,
    2, 0, 0, 1,
    3, 2, 1, -1,
    3, 1, 0, 1,
    3, 0, 0, 1,
    4, 2, 2, -1,
    4, 1, 0, 1,
    4, 0, 0, 1],[],4);
s(end).norm = "N";

% ## ----------
%
% ## numbers follow as production rules in fig 2.16

s(end+1).pattern = 'typeIVc';
s(end).mx = reshape([
    1, 1, 1, -1,
    1, 0, 0, 1,
    2, 1, 2, -1,
    2, 0, 0, 1,
    3, 1, 3, -1,
    3, 0, 0, 1,
    4, 2, 1, -1,
    4, 1, 0, 1,
    4, 0, 0, 1,
    5, 2, 2, -1,
    5, 1, 0, 1,
    5, 0, 0, 1,
    6, 2, 3, -1,
    6, 1, 0, 1,
    6, 0, 0, 1,
    7, 3, 1, -1,
    7, 2, 0, 1,
    7, 1, 0, 1,
    7, 0, 0, 1,
    8, 3, 2, -1,
    8, 2, 0, 1,
    8, 1, 0, 1,
    8, 0, 0, 1,
    9, 3, 3, -1,
    9, 2, 0, 1,
    9, 1, 0, 1,
    9, 0, 0, 1],[],4);
s(end).norm = "N";

% #############################
% ##
% ## Mori's asymmetric step-constrained pattern. Normalized in the
% ## reference length.
% ##
% ## Mori, A.; Uchida, S.; Kurazume, R.; Taniguchi, R.; Hasegawa, T. &
% ## Sakoe, H. Early Recognition and Prediction of Gestures Proc. 18th
% ## International Conference on Pattern Recognition ICPR 2006, 2006, 3,
% ## 560-563
% ##

s(end+1).pattern = 'mori2006';
s(end).mx = reshape([
    1, 2, 1, -1,
    1, 1, 0, 2,
    1, 0, 0, 1,
    2, 1, 1, -1,
    2, 0, 0, 3,
    3, 1, 2, -1,
    3, 0, 1, 3,
    3, 0, 0, 3],[],4);
s(end).norm = "M";

% ## Completely unflexible: fixed slope 1. Only makes sense with
% ## open.begin and open.end
s(end+1).pattern = 'rigid';
s(end).mx = [
    1, 1, 1, -1;
    1, 0, 0, 1];
s(end).norm = "N";

%%
%generates patterns with local slope constraints:
%[1/2, 2], [1/3, 3], ..., [1/N, N]
N = 21;
for slope=2:N
    s(end+1) = gen_symmetric_pattern(slope);
end


%% handle unspecified norms:
for i=1:length(s)
    if isempty(s(i).norm)
        s(i).norm = '';
    end
end

%%
s = struct2table(s);


end

function [s] = gen_symmetric_pattern(slope)
%designed to allow very high or low local slope but not 0/inf
% ## First column: enumerates step patterns.
% ## Second   	 step in query index
% ## Third	 step in reference index
% ## Fourth	 weight if positive, or -1 if starting point

%modeled after symmetricP05:
% s(end).mx = [...
%     1, 1, 3, -1;
%     1, 0, 2, 2;
%     1, 0, 1, 1;
%     1, 0, 0, 1;
%     2, 1, 2, -1;
%     2, 0, 1, 2;
%     2, 0, 0, 1;
%     3, 1, 1, -1;
%     3, 0, 0, 2;
%     4, 2, 1, -1;
%     4, 1, 0, 2;
%     4, 0, 0, 1;
%     5, 3, 1, -1;
%     5, 2, 0, 2;
%     5, 1, 0, 1;
%     5, 0, 0, 1];

s.pattern = ['symmetricS' num2str(slope)];

%define steps algorithmically
mx = [];
c=1;
for i=slope:-1:1
    mx = [mx; c 1 i -1]; %#ok<*AGROW> 
    for j=(i-1):-1:0
        w=1;
        if j==(i-1), w=2; end %diagonals            
        mx = [mx; c 0 j w];        
    end
    c=c+1;
end

%symmetric version by exchanging reference and query dimensions
mxs = [mx(:,1)+max(mx(:,1)) mx(:,[3 2 4])];

mx = [mx; mxs]; %combine
mx = mx(1:end-2,:); %remove duplicate last pattern

s.mx = mx;
s.norm = "N+M";
end