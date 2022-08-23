fprintf(mfilename)
addpath('..');
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
%
M_GC = 2000;
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
addpath(fullfile('..','Code_Parameter_Generation'));
addpath(fullfile('..','Code_Result_Analysis'));
% Load prime result, no Ep masking
p = base_case_parameters();
p.masking = 1; p.steric = 0.3;
result = cell(1,5);

p.vaxnum = 2; p.tmax = 42;
[result{1},~] = load_result(p);

p.vaxnum = 3; p.tmax = 28;
[result{2},~] = load_result(p);

p.earlybooster = 42;
[result{3},~] = load_result(p);

p.earlybooster = 0; vaxnum = 4;
[result{4},~] = load_result(p);

n=4;
summary = cell(1,n);
WT_titer = zeros(2,n);
Var_titer = zeros(2,n);
for i=1:n
    summary{i} = AbConcentration(result{i}, result{i}.param);
    WT_titer(:,i) = (summary{i}.titer_geomean{1}(:,end));
    Var_titer(:,i) = (summary{i}.titer_geomean{2}(:,end));
end
x=1:n;
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'x','WT_titer', 'Var_titer')


