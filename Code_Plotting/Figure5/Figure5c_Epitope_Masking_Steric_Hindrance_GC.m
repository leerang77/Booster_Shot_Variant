fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
result = cell(1,2);

p = base_case_parameters();
p.masking = 1; p.steric = 0.3;
[result{1},~] = load_result(p);
p.vaxnum = 2; p.tmax = 180;
[result{2},~] = load_result(p);
param{1} = result{1}.param;
param{2} = result{2}.param;
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
