fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
[result, param] = loadFigure2Data();
result{2} = result{3}; param{2} = param{3};
for vax=1:2
    [pc{vax},mem{vax}] = MemPCnum(param{vax}, result{vax});
end
time = [param{1}.tspan_summary, param{2}.tspan_summary(1:end)+28];
tidx1 = 1:length(param{1}.tspan_summary);
tidx2 = length(param{1}.tspan_summary)+1:length(time);
for ep=1:2
pcnum{ep} = [pc{1}{ep}, pc{2}{ep}(:,1:end)];
memnum{ep} = [mem{1}{ep}, mem{2}{ep}(:,1:end)];
end
clear result
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile)