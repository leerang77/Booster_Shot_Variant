fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
% Load prime result, no Ep masking
p = base_case_parameters();
p.masking = 1;
p.vaxnum = 2; p.tmax = 180;
sterics = [0, 0.2, 0.3, 0.4, 0.5, 0.7];
n = length(sterics);
result = cell(1,n);
param = cell(1,n);
memnum = zeros(2,n);
for i=1:length(result)
    p.steric = sterics(i);
    [result{i},~] = load_result(p);
    param{i} = result{i}.param;
    buildup3D = @(x) reshape(x, size(x,1), [], param{i}.n_ep);
    memnum_GC = buildup3D(squeeze(result{i}.output.memnumbytime(1:2:end-1, :, 2)));
for ep=1:2
   memnum(ep,i) = mean(memnum_GC(:,end,ep));   
end
result{i} = [];
end
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile)
