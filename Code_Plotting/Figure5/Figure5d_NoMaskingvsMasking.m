fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
result = cell(1,2);
p = base_case_parameters();
p.vaxnum = 3; p.tmax = 28;
[result{1},~] = load_result(p);
p.masking = 1; p.steric = 0.3;
[result{2},~] = load_result(p);
n=2;
summary = cell(1,n);
vax3.WT_titer = zeros(2,n);
vax3.Var_titer = zeros(2,n);
for i=1:n
    summary{i} = AbConcentration(result{i}, result{i}.param);
    vax3.WT_titer(:,i) = (summary{i}.titer_geomean{1}(:,end));
    vax3.Var_titer(:,i) = (summary{i}.titer_geomean{2}(:,end));
end
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'vax3')

