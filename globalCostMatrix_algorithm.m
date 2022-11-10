function [gcm,updates] = globalCostMatrix_algorithm(lm,step_pattern)

nsteps = size(step_pattern.mx,1);

pn = step_pattern.mx(:,1); %pattern id
di = step_pattern.mx(:,2); %delta i
dj = step_pattern.mx(:,3); %delta j
sc = step_pattern.mx(:,4); %step cost

%cost matrix
cm = nan(size(lm));
cm(1,1) = lm(1,1);

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

        updates(i,j).cm = cm;
        updates(i,j).dm = dm;
        updates(i,j).row = i;
        updates(i,j).col = j;
        updates(i,j).costs = costs;
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

