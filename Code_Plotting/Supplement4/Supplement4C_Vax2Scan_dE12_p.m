fprintf(mfilename)
addpath('..');
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
%% Bolus Vax2 Scan - Varying dE12 and p
param = base_case_parameters();
param.vaxnum = 2; param.tmax = 180;
x = [0, 0.2, 0.4, 0.6, 0.8, 1];
y = [0.5,0.25,0.1,0.02];
memnum = repmat({zeros(length(x), length(y))},1,2);
for i=1:length(x)
    for j=1:length(y)
        param.dE12 = x(i);
        param.p = y(j);
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
