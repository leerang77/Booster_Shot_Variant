fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
% Load prime result, no Ep masking
p = cell(1,3);
p{1} = base_case_parameters();
p{2} = p{1}; p{3} = p{1};
p{2}.vaxnum = 2; p{2}.tmax = 180;
p{3}.vaxnum = 3; p{3}.tmax = 28;

result = cell(3,3);
outputprob = [0.05, 0.10, 0.05]; outputpcfrac = [0.1, 0.1, 0.5];
summary = cell(3,3); WT_titer = cell(1,3); Variant_titer = cell(1,3);
tspan = cell(3,3);
titers = {zeros(3,3),zeros(3,3)};
pc = cell(3,3); mem = cell(3,3); memnum = cell(3,2);
for i=1:3
    for j=1:3
        i
        j
        p{j}.outputprob = outputprob(i); p{j}.outputpcfrac = outputpcfrac(i);
        [result{i,j},~] = load_result(p{j});
        summary{i,j} = AbConcentration(result{i,j}, result{i,j}.param);
        tspan{i,j} = result{i,j}.param.tspan_summary;
        [pc{i,j},mem{i,j}] = MemPCnum(result{i,j}.param, result{i,j});
        result{i,j} = [];
    end
WT_titer{i} = sum([summary{i,1}.titer_geomean{1}, summary{i,2}.titer_geomean{1}(:,2:end)]);
Variant_titer{i} = sum([summary{i,1}.titer_geomean{2}, summary{i,2}.titer_geomean{2}(:,2:end)]);

for ep=1:2
memnum{i,ep} = [mem{i,1}{ep}, mem{i,2}{ep}(:,2:end)];
end

for ep=1:2
titers{ep}(i,1) = sum(summary{i,2}.titer_geomean{ep}(:,42*4+1));
titers{ep}(i,2) = sum(summary{i,2}.titer_geomean{ep}(:,150*4+1));
titers{ep}(i,3) = sum(summary{i,3}.titer_geomean{ep}(:,end));
end
end
time = [tspan{1,1}, tspan{1,2}(2:end)+28];



outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'time','WT_titer','Variant_titer','memnum','titers')
