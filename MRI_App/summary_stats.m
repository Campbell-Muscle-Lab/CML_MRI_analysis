function out = summary_stats(x,varargin);
% Returns summary stats

% Defaults
params.ignore_NaNs = 1;

% Update
params = parse_pv_pairs(params,varargin);

% Screen x if required
if (params.ignore_NaNs)
    x = x(~isnan(x));
end

% Error checking
[r,c]=size(x);
if ((r>1)&(c>1))
%     error('summary_stats.m requires a vector');
end

out.n = numel(x);
out.mean = mean(x);
out.sd = std(x);
out.sem = std(x)/sqrt(out.n);
[out.min, out.min_index] = min(x);
[out.max, out.max_index] = max(x);
out.range = max(x) - min(x);