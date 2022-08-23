fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
p = base_case_parameters();
w2s = [0.3, 0.5, 0.7, 1];

memnum = {zeros(2,length(w2s)), zeros(2,length(w2s))};

for idx = 1:length(w2s)
    p.w2 = w2s(idx);
    p.vaxnum = 1; p.tmax = 28;
    [result,~] = load_result(p);
    n = size(result.output.memnumbytime,2)/result.param.n_ep;
    for ep=1:2
        memnum{1}(ep,idx) = mean(squeeze(result.output.memnumbytime(1:2:end-1,n*ep,2)));
    end
    
    p.vaxnum = 2; p.tmax = 180;
    [result,~] = load_result(p);
    n = size(result.output.memnumbytime,2)/result.param.n_ep;
    for ep=1:2
        memnum{2}(ep,idx) = mean(squeeze(result.output.memnumbytime(1:2:end-1,n*ep,2)));
    end
    clear result summary
end
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'w2s', 'memnum')
