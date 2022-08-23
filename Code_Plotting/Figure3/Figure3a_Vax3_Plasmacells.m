fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%% Load prime result, no Ep masking
p = base_case_parameters();
p.vaxnum = 3;
p.tmax = 28;
result = cell(1,1);
affcnts_gcpc = cell(1,1);
affcnts_egcpc = cell(1,1);
[result{1},~] = load_result(p);
param{1} = result{1}.param;
for vaxcase=1:1
    result{vaxcase}.output.finalpc = reshape(permute(result{vaxcase}.output.finalpc,[3,2,1]),param{vaxcase}.pcfieldnum,[]);
    result{vaxcase}.output.finalpc = alignLeft(result{vaxcase}.output.finalpc, 1);
    pc_cumul = [result{vaxcase}.output.finalpc, alignLeft(result{vaxcase}.dead.plasmaCells,1)];
    [~, ctrs, affcnts_gcpc{vaxcase}] = plotHist(pc_cumul);
    if param{vaxcase}.vaxnum > 1
        EGCpc_cumul = [result{vaxcase}.plasmaCellsEGC, alignLeft(result{vaxcase}.dead.plasmaCellsEGC,1)];
        [~, ctrs, affcnts_egcpc{vaxcase}] = plotHist(EGCpc_cumul);
    else
        affcnts_egcpc{vaxcase} = repmat({zeros(1,length(ctrs))},2,2);
    end
end
clear result
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'ctrs', 'affcnts_gcpc', 'affcnts_egcpc')

%%
function [aff, ctrs, affcnts] = plotHist(cellarr)
targetidx = cell(1,2);
    targetidx{1} = cellarr(2,:)==1;
    targetidx{2} = cellarr(2,:)==max(cellarr(2,:));
aff = cell(2,2);
E{1} = cellarr(3,:);
E{2} = cellarr(4,:);
bins = 6.75:0.5:11.25;
ctrs = bins(1:end-1)+0.5*(bins(2)-bins(1));
affcnts = cell(2,2);
for strain=1:2
    for ep=1:2
        idx = targetidx{ep};
        aff{strain,ep}=E{strain}(idx);
        aff{strain,ep}(aff{strain,ep}<4) = 4;
        affcnts{strain,ep} = histcounts(aff{strain,ep},bins);
    end
end
end
