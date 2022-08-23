fprintf(mfilename)
addpath('..');
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
%% Bolus Vax1 Scan - Varying E1h and dE12
% Load prime result, no Ep masking
param = base_case_parameters();
x = 7:0.2:8;
y = 0:0.2:1;
memnum = repmat({zeros(length(x), length(y))},1,2);
for i=1:length(x)
    for j=1:length(y)
        param.E1h = x(i);
        param.dE12 = y(j);
        inputs = {param.vaxnum, param.E1h, param.dE12, param.p, ...
            param.masking, param.C0, param.w1, param.w2, param.steric,...
            param.outputprob, param.outputpcfrac, param.rho, param.tmax};
        inputs = cellfun(@(x) num2str(x), inputs, 'UniformOutput',false);
        fnm = fullfile('..','ParameterScanSummary', [strjoin(inputs(1:end),'_'), '.mat']);
        result = load(fnm);
        for ep=1:2
            memnum{ep}(i,j) = result.memnum{ep};
        end
    end
end
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile)