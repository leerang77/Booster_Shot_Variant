fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
[result, param] = loadFigure2Data();
result{2} = result{3}; param{2} = param{3};
rng(100);
M_GC = 2000;
nn = randsample(M_GC, 100)';
result1 = cell(1,2);
for i=1:2
    result1{i}.gc.numbytime = result{i}.gc.numbytime;
end
result = result1; clear result1
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile)