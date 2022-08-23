fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
% Load Vax3 result, no Ep masking
p = base_case_parameters();
p.vaxnum = 3;
p.tmax = 28;
result = cell(1,1);
[result{1},~] = load_result(p);
n=1;
summary = cell(1,n);
vax3.WT_titer = zeros(2,n);
vax3.Var_titer = zeros(2,n);
for i=1:n
    summary{i} = AbConcentration(result{i}, result{i}.param);
    vax3.WT_titer(:,i) = (summary{i}.titer_geomean{1}(:,end));
    vax3.Var_titer(:,i) = (summary{i}.titer_geomean{2}(:,end));
end
% Load Vax2 result, no Ep masking
p.vaxnum = 2;
p.tmax = 180;
result = cell(1,1);
[result{1},~] = load_result(p);
n=1;
summary = cell(1,n);
vax2.WT_titer = zeros(2,n);
vax2.Var_titer = zeros(2,n);
for i=1:n
    summary{i} = AbConcentration(result{i}, result{i}.param);
    vax2.WT_titer(:,i) = (summary{i}.titer_geomean{1}(:,42*4+1));
    vax2.Var_titer(:,i) = (summary{i}.titer_geomean{2}(:,42*4+1));
    vax2long.WT_titer(:,i) = (summary{i}.titer_geomean{1}(:,150*4+1));
    vax2long.Var_titer(:,i) = (summary{i}.titer_geomean{2}(:,150*4+1));
end
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'vax2', 'vax2long', 'vax3')
