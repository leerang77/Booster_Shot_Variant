fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
result = cell(1,2);

p = base_case_parameters();
p.masking = 1; p.steric = 0.3;
p.vaxnum = 2; p.tmax = 180;
[result{2},~] = load_result(p);

param{2} = result{2}.param;

conc = cellfun(@(obj) obj.concarray_Epmask, result{2}.conc, 'UniformOutput', false);
conc = vertcat(conc{:});
conc = {squeeze(geomean(conc(1:2:end-1,:,:),1)); squeeze(geomean(conc(2:2:end,:,:),1))};
Ceff = (param{2}.Ageff*conc{1} + conc{2});

outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'Ceff', 'param')
