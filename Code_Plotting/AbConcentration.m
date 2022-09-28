function summary = AbConcentration(result, param)
%% Documentation
% Summary: Obtain summary statistics of antibody titers. This function is
%   called by several plotting functions. 
% Inputs:
%   result - struct containing the simulation results (all repeats)
%   param - struct containing the simulation parameters
% Output:
%   summary - struct containing the summary statistics of antibody titers

n_sim = length(result.conc);
xlength = length(param.tspan_summary);
KaType = {'Kaarray', 'Kaarray_var'};
titers = repmat({zeros(param.n_ep*n_sim, xlength)},1,2);
[titer_mean,titer_std,titer_geomean,titer_geostd,titerfraction] = ...
    deal(cell(1,2));
affinity_weighted = 1;

for variant=1:2
    for i=1:n_sim
        if affinity_weighted
            titer = (result.conc{i}.concarray(2:end,2:end,:).*...
                result.conc{i}.(KaType{variant}));
        else
            titer = result.conc{i}.concarray(2:end,2:end,:);
        end
        titer = squeeze(sum(titer,1));
        titers{variant}(param.n_ep*(i-1)+1:param.n_ep*i,:) = titer;
    end
    titer_mean{variant} = ...
        squeeze(mean(reshape(titers{variant},param.n_ep,n_sim,[]),2));
    titer_std{variant} = ...
        squeeze(std(reshape(titers{variant},param.n_ep,n_sim,[]),0,2));
    titer_geomean{variant} = ...
        squeeze(geomean(reshape(titers{variant},param.n_ep,n_sim,[]),2));
    titer_geostd{variant} = ...
        squeeze(geostd(reshape(titers{variant},param.n_ep,n_sim,[]),2));
    titerfraction{variant} = ...
        titer_mean{variant}(param.n_ep,end)/titer_mean{variant}(1,end);

    summary.titer_mean = titer_mean;
    summary.titerfraction = titerfraction;
    summary.titer_std = titer_std;
    summary.titer_geomean = titer_geomean;
    summary.titer_geostd = titer_geostd;
end
end

%% Subfunction
function geosig = geostd(X,dim)
% Summary: Obtains the geometric standard deviation of X along the given
%          dimension.
    Y = log(X);
    n = size(X,dim);
    mu = mean(Y,dim);
    geosig = exp(sqrt(sum((Y-mu).^2,dim)/(n-1)));
end