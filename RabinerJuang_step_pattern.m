function [s] = RabinerJuang_step_pattern(ptype,slope_weighting,smoothed)

if nargin==0
    ptype = 1;
    slope_weighting = 'a';
    smoothed = false;
end

switch(ptype)
    case 1
        steps = {
            [1 0];
            [1 1];
            [0 1]};

    case 2
        steps = {
            [1 1; 1 0];
            [1 1];
            [1 1; 0 1]};        

    case 3
        steps = {
            [2 1];
            [1 1];
            [1 2]};          

    case 4
        steps = {
            [1 1; 1 0];
            [1 2; 1 0];
            [1 1];
            [1 2]};          

    case 5
        steps = {
            [1 1; 1 0; 1 0];
            [1 1; 1 0];
            [1 1];
            [1 1; 0 1];
            [1 1; 0 1; 0 1]};          

    case 6
        steps = {
            [1 1; 1 0; 1 0];
            [1 1];
            [1 1; 0 1; 0 1]};          

    case 7 
        steps = {
            [1 1; 1 0; 1 0];
            [1 2; 1 0; 1 0];
            [1 3; 1 0; 1 0];
            [1 1; 1 0];
            [1 2; 1 0];
            [1 3; 1 0];
            [1 1];
            [1 2];
            [1 3]};         

end

%%
s.pattern = ['RJType' num2str(ptype) slope_weighting];

mx = [];
for j=1:size(steps,1)
    P(j) = initP(j,slope_weighting,smoothed);
    for k=1:size(steps{j},1)
        P(j).i = [P(j).i steps{j}(k,1)];
        P(j).j = [P(j).j steps{j}(k,2)];
    end
    mx = [mx; getmx(P(j))];
end
s.mx = mx;

s.norm = {};
switch(slope_weighting)
    case 'c'
        s.norm = {'N'};
    case 'd'
        s.norm = {'N+M'};
end

end

%%
function [mx] = getmx(P)

ia = P.i;
ja = P.j;

si = cumsum(ia);
sj = cumsum(ja);

ni = max(si) - si;
nj = max(sj) - sj;

switch(P.subtype)
    case 'a'
        w = min([ia; ja]);
    case 'b'
        w = max([ia; ja]);
    case 'c'
        w = ia;
    case 'd'
        w = ia + ja;
end

if P.smoothing
    w(2:end) = mean(w(2:end));
end

w(1) = -1.0;

mx = zeros(length(P.i),4);
mx(:,1) = P.pid;
mx(:,2) = ni';
mx(:,3) = nj';
mx(:,4) = w';

end


%%
function [P] = initP(pid,subtype,smoothing)
P.subtype = subtype;
P.smoothing = smoothing;
P.pid = pid;
P.i = [1];
P.j = [1];
end


