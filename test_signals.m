function [s] = test_signals(varargin)

% This function defines a variety of signals used for testing.
%   Not every signal defined below is used.

dbstop if error;

p = inputParser;

def_name = {};
def_params = {};

addOptional(p,'name',def_name,@(x)ischar(x) || iscell(x));
addOptional(p,'params',def_params,@(x)iscell(x));

parse(p,varargin{:});

r = p.Results;
if ischar(r.name)
    r.name = {r.name};
end

params = p.Results.params;

%% signal data function handles
F.lin_rate_var          = @()lin_rate_var(params{:});
F.lin_rate_var2         = @()lin_rate_var2(params{:});
F.many_to_many_example  = @()many_to_many_example(params{:});
F.alignment_example     = @()alignment_example(params{:});
F.alignment_example2    = @()alignment_example2(params{:});
F.example_1             = @()example_1(params{:});
F.cos_freq_var          = @()cos_freq_var(params{:});
F.cos_freq_var_trunc    = @()cos_freq_var_trunc(params{:});
F.cos_freq_pert         = @()cos_freq_pert(params{:});
F.cos_amp_var           = @()cos_amp_var(params{:});
F.cos_per_incom         = @()cos_per_incom(params{:});
F.cos_phase_pert        = @()cos_phase_pert(params{:});
F.growth_decay_constant = @()growth_decay_constant(params{:});
F.growth_decay_variable = @()growth_decay_variable(params{:});
F.growth_decay_set      = @()growth_decay_set(params{:});
F.growth_decay_example  = @()growth_decay_example(params{:});
F.growth_decay_phases   = @()growth_decay_phases(params{:});
F.growth_incom          = @()growth_incom(params{:});
F.growth_set            = @()growth_set(params{:});
F.local_rate_example    = @()local_rate_example(params{:});
F.similar_shape         = @()similar_shape(params{:});
F.slope_constraint_example  = @()slope_constraint_example(params{:});
F.edge_effect_example       = @()edge_effect_example(params{:});
F.global_window_example     = @()global_window_example(params{:});
F.open_endpoint_example     = @()open_endpoint_example(params{:});


%%

ff = fieldnames(F);

for i=1:length(ff)
    S{i} = F.(ff{i})();
    S{i}.name = ff{i};
end

%%
for i=1:length(S)
    T(i).name = S{i}.name;
    T(i).S = S{i};
end

T = struct2table(T);

if nargin==0
    s = T;
else
    s = T.S{ismember(T.name,r.name)};
end

end

%%
function [S] = local_rate_example(varargin)

dt = 0.025;
t = 0:dt:0.5;
S.rate = 2;
S.rate_range = S.rate + 0.5*[1 -1];

rates{1} = S.rate*ones(size(t));
rates{2} = linspace(S.rate_range(1),S.rate_range(2),length(t));

for i=1:2
    S.dt(i) = dt;
    S.t{i} = t;
    S.rates{i} = rates{i}*dt;
    S.state{i} = [0 cumsum(S.rates{i}(1:end-1))];
    S.X{i} = exp(-(2)*S.state{i});
end

end


%%
function [S] = growth_decay_constant(varargin)
x0 = 0;
sigma0 = 1;
dt = 0.001;
S.dt = dt;
S.constant_rates = [0.5 1 1.5];

for i=1:length(S.constant_rates)
    S.rates{i} = S.constant_rates(i)*ones(1,10000)*dt;
    S.X{i} = gd_process(x0,sigma0,S.rates{i});
    S.t{i} = 0:dt:(length(S.X{i})-1)*dt;

    S.rates{i} = S.rates{i}(1:length(S.t{i}));
end

end

%%
function [S] = growth_decay_variable(varargin)
x0 = 0;
sigma0 = 1;
dt = 0.001;

S.rate_var_fac = [1 0 -2];
S.dt = dt;

for i=1:length(S.rate_var_fac)

    c = 1;
    m = S.rate_var_fac(i);

    S.rates{i} = linspace(c-m,c+m,2000)*dt;

    S.X{i} = gd_process(x0,sigma0,S.rates{i});
    S.t{i} = 0:dt:(length(S.X{i})-1)*dt;
    S.rates{i} = S.rates{i}(1:length(S.t{i}));
end

end

%%
function [S] = growth_decay_example(varargin)
x0 = 0;
sigma0 = 1;
dt = 0.001;

S.rate_var_fac = [1 0 -1];
S.dt = dt;

for i=1:length(S.rate_var_fac)

    c = 2;
    m = S.rate_var_fac(i);

    S.rates{i} = linspace(c-m,c+m,1200)*dt;

    S.X{i} = gd_process(x0,sigma0,S.rates{i});
    S.t{i} = 0:dt:(length(S.X{i})-1)*dt;
    S.rates{i} = S.rates{i}(1:length(S.t{i}));
end

end

%%
function [S] = growth_set(varargin)

p = inputParser;
def_n = 5;
addParameter(p,'n',def_n);
parse(p,varargin{:});

n = p.Results.n;

x0 = 0;
dt = 0.001;

%S.constant_rates = exp(linspace(log(S.rate_range(1)),log(S.rate_range(2)),n));
S.constant_rates = 2*(1/2).^linspace(1,-1,n);
S.dt = dt;

for i=1:length(S.constant_rates)
    [S.X{i},S.rates{i}] = growth_process(x0,S.constant_rates(i)*dt);
    S.t{i} = 0:dt:(length(S.X{i})-1)*dt;
end

end

%%
function [x,rate] = growth_process(x0,rate)
x = x0;
i = 1;
if numel(rate)==1, rate = rate*ones(1,10000); end
while x(i)<=1
    i=i+1;
    x(i) = x(i-1)+rate(i);
end
x = x(x<=1);
rate = rate(1:length(x));
end

%%
function [S] = growth_decay_set(varargin)

p = inputParser;
def_n = 5;
addParameter(p,'n',def_n);
parse(p,varargin{:});

n = p.Results.n;

x0 = 0;
sigma0 = 1;
dt = 0.001;

S.rate_range = [1/3 3];
S.phase_rates = exp(linspace(log(S.rate_range(1)),log(S.rate_range(2)),n));
S.phase_rates = [S.phase_rates; fliplr(S.phase_rates)]';
S.dt = dt;

for i=1:size(S.phase_rates,1)
    [S.X{i},S.rates{i}] = gd_process_phase(x0,sigma0,S.phase_rates(i,:)*dt);
    S.t{i} = 0:dt:(length(S.X{i})-1)*dt;
    S.rates{i} = S.rates{i}(1:length(S.t{i}));
end

end

%%
function [S] = growth_decay_phases(varargin)

p = inputParser;
def_n = 5;
addParameter(p,'n',def_n);
parse(p,varargin{:});

n = p.Results.n;

x0 = 0;
sigma0 = 1;
dt = 0.001;

S.rate_range = [1*[1 3/2]; 1*[1 2]];
S.phase_rates = [
    exp(linspace(log(S.rate_range(1,1)),log(S.rate_range(1,2)),n));
    exp(linspace(log(S.rate_range(2,1)),log(S.rate_range(2,2)),n))]';
S.dt = dt;

for i=1:size(S.phase_rates,1)
    [S.X{i},S.rates{i}] = gd_process_phase(x0,sigma0,S.phase_rates(i,:)*dt);
    S.t{i} = 0:dt:(length(S.X{i})-1)*dt;
    S.rates{i} = S.rates{i}(1:length(S.t{i}));
end

end

%%
function [x] = gd_process(x0,sigma0,rates)
x = x0;
sigma = sigma0;
for i=2:length(rates)
    x(i) = x(i-1)+sigma*rates(i);

    if x(i)>=1, sigma = -1; end
    if x(i)<=0, break; end
end
x = x(x>=0);
end

%%
function [x,rate] = gd_process_phase(x0,sigma0,rates)
x = x0;
sigma = sigma0;
rate = rates(1);
i=1;
while 1
    i=i+1;
    rate(i) = rate(i-1);    
    x(i) = x(i-1)+sigma*rate(i-1);

    if x(i)>=1, sigma = -1; rate(i) = rates(2); end
    if x(i)<=0, break; end
end
x = x(x>=0);
end

%% cosine (1-period) frequency variation
function [S] = cos_freq_var(varargin)

p = inputParser;
def_n = 5;
def_n_periods = 1;
addParameter(p,'n',def_n);
addParameter(p,'n_periods',def_n_periods);
parse(p,varargin{:});

n = p.Results.n;
n_per = p.Results.n_periods;

st = dbstack;
S.name = st.name;
S.descr = 'cosines of varying frequency';
S.frq = 1./(linspace(0.2,1,n));
S.dt = 1e-3;

for i=1:length(S.frq)
    S.t{i} = (0:S.dt:(n_per/S.frq(i)));
    S.X{i} = cos(2*pi*S.frq(i)*S.t{i});
end

S.formula = '$$y=\cos(2\pi f_i t)$$';
S.labels = arrayfun(@(c,d){sprintf('f_{%i}=%1.2f',c,d)},1:length(S.frq),S.frq);

end

%% cosine (1-period) frequency variation, truncated
function [S] = cos_freq_var_trunc(varargin)

p = inputParser;
def_n = 5;
addParameter(p,'n',def_n);
parse(p);
n = p.Results.n;

st = dbstack;
S.name = st.name;
S.descr = 'cosines of varying frequency, truncated';
S.frq = 1./(linspace(0.1,1,n));
S.dt = 1e-3;

for i=1:length(S.frq)
    S.t{i} = (0:S.dt:(1/S.frq(i)));
    S.X{i} = cos(2*pi*S.frq(i)*S.t{i});
end

%truncate
max_len = length(S.X{ceil(length(S.frq)/2)});
for i=1:length(S.frq)
    ixs = 1:min(max_len,length(S.X{i}));
    S.t{i} = S.t{i}(ixs);
    S.X{i} = S.X{i}(ixs);
end

S.formula = '$$y=\cos(2\pi f_i t)$$';
S.labels = arrayfun(@(c,d){sprintf('f_{%i}=%1.2f',c,d)},1:length(S.frq),S.frq);

end

%% linear (1-period) growth rate variation
function [S] = lin_rate_var(varargin)

p = inputParser;
def_n = 5;
addParameter(p,'n',def_n);
parse(p,varargin{:});
n = p.Results.n;

st = dbstack;
S.name = st.name;
S.descr = 'linear growth-decay at constant rate';
S.frq = 1./(linspace(0.5,1,n));
S.dt = 1e-3;

for i=1:length(S.frq)
    S.t{i} = (0:S.dt:(1/S.frq(i)));
    S.X{i} = 1 - 2*S.frq(i)* abs(S.t{i}-(1/(2*S.frq(i))));
end

S.formula = '$$y=1-2r\bigg|t-\frac{1}{2r}\bigg|$$';
S.labels = arrayfun(@(c,d){sprintf('r_{%i}=%1.2f',c,d)},1:length(S.frq),S.frq);

end

%% linear (1-period) growth rate variation
function [S] = lin_rate_var2(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'linear growth-decay at constant rate';
S.frq = [0.5 1 2.0];
S.dt = 1e-3;

for i=1:length(S.frq)
    S.t{i} = (0:S.dt:(1/S.frq(i)));
    S.X{i} = 1 - 2*S.frq(i)* abs(S.t{i}-(1/(2*S.frq(i))));
end

S.formula = '$$y=1-2r\bigg|t-\frac{1}{2r}\bigg|$$';
S.labels = arrayfun(@(c,d){sprintf('r_{%i}=%1.2f',c,d)},1:length(S.frq),S.frq);

end

%% example of many-to-many mapping
function [S] = many_to_many_example(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'signals with many-to-many mapping';
S.dt = 1;

S.t{1} = 1:8;
S.X{1} = [0 1 2 3 4 4 4 5]/5;
S.t{2} = 1:8;
S.X{2} = [0 1 1 1 2 3 4 5]/5;

S.formula = '';
S.labels = {'s1','s2'};

end

%% example of many-to-many mapping
function [S] = alignment_example(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'signals with many-to-many mapping';
S.dt = 1;

S.t{1} = 1:8;
S.X{1} = [0 1 2 3 3.90 4 4.1 5]/5;
S.t{2} = 1:8;
S.X{2} = [0 0.90 1 1.1 2 3 4 5]/5;

S.formula = '';
S.labels = {'s1','s2'};

end

%% example of many-to-many mapping
function [S] = alignment_example2(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'example for illustration of distance';
S.dt = 1;

S.t{1} = 1:8;
S.X{1} = [0 1 2 3 4 4 5 6]/6;
S.t{2} = 1:9;
S.X{2} = [0 1 1 1 2 3 4 5 5]/5;

S.formula = '';
S.labels = {'s1','s2'};

end
%% cosine frequency variation (only three)
function [S] = cos_freq_var_simple(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'cosines of varying frequency';
S.frq = 1./([0.5 0.75 1]);
S.dt = 1e-3;

for i=1:length(S.frq)
    S.t{i} = (0:S.dt:(1/S.frq(i)));
    S.X{i} = cos(2*pi*S.frq(i)*S.t{i});
end

S.formula = '$$y=\cos(2\pi f_i t)$$';
S.labels = arrayfun(@(c,d){sprintf('f_{%i}=%1.0f',c,d)},1:length(S.frq),S.frq);

end

%% cosine frequency variation (only three)
function [S] = example_1(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'example signals';
S.frq = 1./([0.5 1]);
S.dt = 1e-3;

for i=1:length(S.frq)
    S.t{i} = (0:S.dt:(1/S.frq(i)));
    S.X{i} = cos(2*pi*S.frq(i)*S.t{i});
end

S.formula = '$$y=\cos(2\pi f_i t)$$';
S.labels = arrayfun(@(c,d){sprintf('f_{%i}=%1.0f',c,d)},1:length(S.frq),S.frq);

end

%% cosine (1-period) phase shifts
function [S] = cos_phase_pert(varargin)

p = inputParser;
def_n = 5;
addParameter(p,'n',def_n);
parse(p);
n = p.Results.n;

st = dbstack;
S.name = st.name;
S.descr = 'cosines with phase offsets';
S.frq = 5;
S.dt = 1e-3;
S.phase_pert = linspace(0,pi,n);

t = (0:S.dt:(1/S.frq - S.dt));

for i=1:length(S.phase_pert)
    S.t{i} = t;
    S.X{i} = cos(2*pi*S.frq.*S.t{i} + S.phase_pert(i));
end

S.formula = '$$\cos(2\pi (f+\sin(2\pi t+\phi_i)) t)$$';
%S.labels = arrayfun(@(c){sprintf('\\omega=%1.0f',c)},S.frq);

end

%% cosine (1-period) frequency perturbations
function [S] = cos_freq_pert(varargin)

p = inputParser;
def_n = 5;
addParameter(p,'n',def_n);
parse(p);
n = p.Results.n;

st = dbstack;
S.name = st.name;
S.descr = 'cosines with perturbed frequencies';
S.frq = n;
S.dt = 1e-3;

%center time of perturbation
S.cycle_pert = 1:n;

sigma = 0.05; %standard deviation
gaussfcn = @(mu)exp(-((t-mu).^2)/(2*sigma^2));

for i=1:length(S.cycle_pert)
    fp = S.frq + 1;
    t = (0:S.dt:(1/S.frq - S.dt));
    tp = (0:S.dt:(1/fp - S.dt));
    x = cos(2*pi*S.frq*t);
    xp = cos(2*pi*fp*tp);
    S.X{i} = [repmat(x,1,S.cycle_pert(i)-1) xp repmat(x,1,n-S.cycle_pert(i))];
    S.t{i} = (0:length(S.X{i})-1)*S.dt;
end

%S.formula = '$$\cos(2\pi (f+a_i \sin(2\pi(2f) t)) t)$$';
S.formula = '$$c_i: \mathrm{cycle\ of\ perturbation}$$';
S.labels = arrayfun(@(c,d){sprintf('c_{%i}=%1.0f',c,d)},1:length(S.cycle_pert),S.cycle_pert);

end

%% cosines with varying amplitude 
function [S] = cos_amp_var(varargin)

p = inputParser;
def_n = 5;
addParameter(p,'n',def_n);
parse(p);
n = p.Results.n;

st = dbstack;
S.name = st.name;
S.descr = 'cosines with amplitude variation';
S.frq = 1;
S.dt = 1e-3;
S.amp_pert = linspace(1,2,n);

t = (0:S.dt:(1/S.frq));

for i=1:length(S.amp_pert)
    S.t{i} = t;
    S.X{i} = S.amp_pert(i)*cos(2*pi*S.frq*S.t{i});
end

S.formula = '$$y=a_i \cos(2\pi f t)$$';
S.labels = arrayfun(@(c,d){sprintf('a_{%i}=%1.1f',c,d)},1:length(S.amp_pert),S.amp_pert);

end

%% cosines with varying number of periods 
function [S] = cos_per_incom(varargin)

p = inputParser;
def_n = 3;
addParameter(p,'n',def_n);
parse(p);
n = p.Results.n;

st = dbstack;
S.name = st.name;
S.descr = 'cosines with amplitude variation';
S.frq = [0.9 0.75 1];
S.dt = 1e-3;
S.periods = [2 2 2.75];
S.phase_offsets = [-1 0 1]*0.01;

for i=1:length(S.periods)
    S.t{i} = (0:S.dt:(S.periods(i)/S.frq(i)));
    S.X{i} = hamming(length(S.t{i}))'.*cos(2*pi*S.frq(i)*(S.t{i}+S.phase_offsets(i)));
end

% S.formula = '$$y=\cos(2\pi f t)$$';
% S.labels = arrayfun(@(c,d){sprintf('a_{%i}=%1.1f',c,d)},1:length(S.amp_pert),S.amp_pert);

end

%% cosines with varying number of periods 
function [S] = growth_incom(varargin)

p = inputParser;
def_n = 3;
addParameter(p,'n',def_n);
parse(p);
n = p.Results.n;

st = dbstack;
S.name = st.name;
S.descr = 'cosines with amplitude variation';
S.frq = 1;
S.dt = 1e-3;
S.rates = [-20 -5 1 1];
S.fcn = {@(t,r)exp(r*t), @(t,r)exp(r*t), @(t,r)1-r*t, @(t,r)exp(r*t)};

for i=1:length(S.rates)
    S.t{i} = (0:S.dt:1);
    S.X{i} = S.fcn{i}(S.t{i},S.rates(i));
end

% S.formula = '$$y=\cos(2\pi f t)$$';
% S.labels = arrayfun(@(c,d){sprintf('a_{%i}=%1.1f',c,d)},1:length(S.amp_pert),S.amp_pert);

end

%% cosines with varying number of periods 
function [S] = similar_shape(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'cosines with transient frequency perturbation';
S.dt = 1e-3;
S.rates = [20 20 20];
S.pert_amps = [0 0.02 0.5];

sigmoidfcn = @(t,mu,b)1./(1+exp(-b*(t-mu)));

sigma = 0.1;
gaussfcn = @(t,mu)exp(-((t-mu).^2)/(2*sigma^2));

for i=1:length(S.rates)
    t = 0:S.dt:1;
    part = sigmoidfcn(t,0.5,S.rates(i));
    S.t{i} = [t t+1];
    S.X{i} = [part fliplr(part)] - S.pert_amps(i)*gaussfcn(S.t{i},1);
end

end

%% example for illustrating use of local slope constraints
function [S] = slope_constraint_example(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'noise that gives bad alignment in standard DTW';
S.dt = 1e-3;
S.frq = [1 1];
S.noise_amp = [0 0.2];
S.noise_loc = 0.50;
S.sigma = 0.05;

gaussfcn = @(t,mu)exp(-((t-mu).^2)/(2*S.sigma^2));

for i=1:length(S.noise_amp)
    S.t{i} = 0:S.dt:1;
    S.X{i} = cos(2*pi*S.frq(i)*S.t{i}) + S.noise_amp(i)*gaussfcn(S.t{i},S.noise_loc);
end

end

%% example for illustrating use of local slope constraints
function [S] = edge_effect_example(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'sample-edge values do not match';
S.dt = 1e-3;
S.begin_val = [0.2 0];
S.end_val = [0.8 1];

for i=1:length(S.begin_val)
    S.t{i} = 0:S.dt:1;
    S.X{i} = linspace(S.begin_val(i),S.end_val(i),length(S.t{i}));
end

end

%% example for illustrating use of local slope constraints
function [S] = global_window_example(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'examples for illustrating global windows';
S.dt = 1e-3;

S.amp_mod = [5 1; 1 5];
S.amp_locs = [0.25 0.75];
S.sigma = 0.025;

gaussfcn = @(t,mu)exp(-((t-mu).^2)/(2*S.sigma^2));

for i=1:size(S.amp_mod,1)
    S.t{i} = 0:S.dt:1;
    S.X{i} = S.amp_mod(i,1)*gaussfcn(S.t{i},S.amp_locs(1)) + ...
             S.amp_mod(i,2)*gaussfcn(S.t{i},S.amp_locs(2));
end

end

%% example for illustrating use of open end/begin
function [S] = open_endpoint_example(varargin)

st = dbstack;
S.name = st.name;
S.descr = 'examples for open endpoint alignments';
S.dt = 1e-2;
S.t_extra = [0 0.2];

for i=1:length(S.t_extra)
    S.t{i} = (0-S.t_extra(i)):S.dt:(1+S.t_extra(i));
    S.X{i} = cos(2*pi*S.t{i});
end

end