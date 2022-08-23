fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%
% Load Vax3 result, no Ep masking
p = base_case_parameters();
p.w2 = 1;
w1s = [0, 500, 100, 20];

memnum = {zeros(2,length(w1s)), zeros(2,length(w1s))};

for idx = 1:length(w1s)
    p.w1 = w1s(idx);
    p.vaxnum = 1; p.tmax = 28;
    [result,~] = load_result(p);
    n = size(result.output.memnumbytime,2)/result.param.n_ep;
    for ep=1:2
        memnum{1}(ep,idx) = mean(squeeze(result.output.memnumbytime(1:2:end-1,n*ep,3)));
    end
    
    p.vaxnum = 2; p.tmax = 180;
    [result,~] = load_result(p);
    n = size(result.output.memnumbytime,2)/result.param.n_ep;
    for ep=1:2
        memnum{2}(ep,idx) = mean(squeeze(result.output.memnumbytime(1:2:end-1,n*ep,3)));
    end
    clear result summary
end
w1s(1) = Inf;
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'w1s', 'memnum')