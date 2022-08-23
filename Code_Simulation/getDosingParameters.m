function [agconc, abconc, Ka, param] = getDosingParameters(param)
%% Documentation
% Summary:
% Initialize the antigen and antibody concentrations, 
% antibody affinity, and the parameters that define the dosing profile

% Outputs:
%   agconc: 1x3 vector; Concentrations of Antigens; 
%                       soluble, IC1, IC2
%   abconc: 3x2 vector; Concentration of antibodies
%           Dim1 - IgM(natural),IgM(immune),IgG; Dim2-Epitopes
%   Ka: 3x2 vector; Antibody binding affinities (Ka)
%   param: parameter struct

% Input:
%   param: parameter struct

%%
% Initialization
[T,k,numshot] = deal(param.T, param.k, param.numshot);
agconc = zeros(1,param.n_ep+1); %
abconc = zeros(3, param.n_ep);
abconc(1,:) = param.IgM0/param.n_ep;
Ka = 10^-3*ones(3,param.n_ep); %nM^-1

if T==0 %bolus
    param.F0 = 0;
    param.dose_t = 0;
    param.dose = param.Ag0;
elseif T>0 %slow delivery. NOT USED IN THIS STUDY
    if numshot == 0
        if k~=0 %increasing dose
            options = optimoptions('fsolve','Display','off');
            param.F0 = fsolve(@(F0) param.Ag0-(F0/param.k*...
                (exp(param.k*param.T)-1)), 0.1, options);
        else %constant dose
            param.F0 = param.Ag0/param.T;
        end   
    elseif numshot > 0
        param.F0 = 0;
        param.dose_t = linspace(0,param.T,numshot+1);
        param.dose_t = param.dose_t(1:end-1);
        param.dose = param.Ag0*exp((0:(numshot-1))*k)...
            /sum(exp((0:(numshot-1))*k));
    end
end
end
