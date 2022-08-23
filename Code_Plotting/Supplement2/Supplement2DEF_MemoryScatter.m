fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
%%
[result, param] = loadFigure2Data();
for vaxcase=1:3
    result{vaxcase}.output.finalmem = reshape(permute(result{vaxcase}.output.finalmem,[3,2,1]),param{vaxcase}.memfieldnum,[]);
    result{vaxcase}.output.finalmem = alignLeft(result{vaxcase}.output.finalmem, 1);
    [memaff{vaxcase}, ctrs, affcnts_gcmem{vaxcase}] = plotHist(result{vaxcase}.output.finalmem);
    if vaxcase > 1
    [memaff_egc{vaxcase}, ctrs, affcnts_egcmem{vaxcase}] = plotHist(result{vaxcase}.finalmemEGC);
    else
        affcnts_egcmem{vaxcase} = repmat({zeros(1,length(ctrs))},2,2);
    end
end
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'memaff','memaff_egc')


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
        aff{strain,ep}(aff{strain,ep}<6) = 6;
        affcnts{strain,ep} = histcounts(aff{strain,ep},bins);
    end
end
end

