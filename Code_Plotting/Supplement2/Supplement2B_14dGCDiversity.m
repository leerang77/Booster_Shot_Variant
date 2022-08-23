fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
%%
[result, param] = loadFigure2Data();
result{2} = result{3}; param{2} = param{3};
M_GC = 2000;
time = 14;
timeidx = (time/7)+1;
clonality = zeros(M_GC,length(time));
bins = linspace(0,1,21);
ctrs = bins(1:end-1)+0.5*(bins(2)-bins(1));
clonality_hist = zeros(length(bins)-1,length(time));
for vaxcase=1:2
    a = squeeze(result{vaxcase}.gc.numbytime(:,113,1));
    gc_success_idx = a>800;
    bcellnum = squeeze(result{vaxcase}.gc.numbylineage(:,timeidx,:));
    clonality(:,vaxcase) = max(bcellnum,[],2)./sum(bcellnum,2);
    clonality_hist(:,vaxcase) = histcounts(clonality(gc_success_idx,vaxcase), bins);
end
clonality_hist = clonality_hist./sum(clonality_hist);
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'ctrs', 'clonality_hist', 'param')
