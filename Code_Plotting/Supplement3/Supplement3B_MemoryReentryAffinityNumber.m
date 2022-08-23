fprintf(mfilename)
addpath('..');
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
% Load prime result, no Ep masking
param = base_case_parameters();
param.vaxnum = 2; param.tmax = 180;
param.memToGCFrac = 0.04;
result = load_result(param);
%% Percentage of memory-derived GC B cells
[memMemoryFraction, maxMemAff, ctrs, vals] = MemoryReactivationStats(result);
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'memMemoryFraction','maxMemAff','ctrs','vals','param')

function [memMemoryFraction, maxMemAff, ctrs, vals] = MemoryReactivationStats(result)
output = result.output.finalmem;
n_GC = size(output,1);
memMemoryFraction = zeros(n_GC,1);
for k=1:n_GC
    lineages = squeeze(output(k,:,1));
    memMemoryFraction(k) = sum(lineages>2010)/sum(lineages>0);
end
maxMemAff = zeros(n_GC,1);
naive = result.naive;
for k=1:n_GC
    maxMemAff(k) = max(naive(k,2100:end,3));
end
bins = 6:0.5:11;
ctrs = bins(1:end-1)+(bins(2)-bins(1))/2;
vals = zeros(size(ctrs));
for i=1:length(vals)
    idx = find(maxMemAff>bins(i) & maxMemAff<bins(i+1));
    vals(i) = mean(memMemoryFraction(idx));
end
end
% figure
% histogram(memMemoryFraction, 0:0.05:1)
% xlabel('Fraction of memory-derived GC B cells')
% ylabel('Number of GCs')
% set(gcf, 'Position', [100, 100, 300, 300]);