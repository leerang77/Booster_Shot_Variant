fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%
[result, param] = loadFigure2Data();
result{2} = result{3}; param{2} = param{3};
for i=1:2
    summary{i} = AbConcentration(result{i}, param{i});
end
time = [param{1}.tspan_summary, param{2}.tspan_summary(2:end)+28];
WT_titer = [summary{1}.titer_geomean{1}, summary{2}.titer_geomean{1}(:,2:end)];
Variant_titer = [summary{1}.titer_geomean{2}, summary{2}.titer_geomean{2}(:,2:end)];
clear result
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile)


