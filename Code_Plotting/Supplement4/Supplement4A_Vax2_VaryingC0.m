fprintf(mfilename)
matFileName = mfilename;
addpath('..');
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
% Load prime result, no Ep masking
param = base_case_parameters();
param.vaxnum = 2; param.tmax = 180;
C0s = param.C0 * [8, 4, 2, 1, 1/2, 1/4, 1/8];
n = length(C0s);
memnum = zeros(2,n);
for i=1:n
param.C0 = C0s(i);
inputs = {param.vaxnum, param.E1h, param.dE12, param.p, ...
 param.masking, param.C0, param.w1, param.w2, param.steric,...
 param.outputprob, param.outputpcfrac, param.rho, param.tmax};
inputs = cellfun(@(x) num2str(x), inputs, 'UniformOutput',false);
fnm = fullfile('..','ParameterScanSummary', [strjoin(inputs(1:end),'_'), '.mat']);
result = load(fnm);
    for ep=1:2
        memnum(ep,i) = result.memnum{ep};
    end
end
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
outputFile = fullfile('..','Results',[matFileName, '.mat']);
save(outputFile)
