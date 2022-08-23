function summarizeMemoryNum(varargin)
addpath(fullfile('..','Code_Parameter_Generation'));
addpath(fullfile('..','Code_Result_Analysis'));
addpath(fullfile('..','Code_Simulation'));

if ~exist('ParameterScanSummary','dir')
    mkdir('ParameterScanSummary')
end

p = initializeParameters(varargin{:});
[result,dirname] = load_result(p);

[~,dirname,~] = fileparts([dirname,'.dummy']);
outputFile = [fullfile('ParameterScanSummary',dirname),'.mat'];
% if exist(outputFile, 'file')
%     fprintf('Output file already exists')
%     return
% end

n = size(result.output.memnumbytime,2)/result.param.n_ep;
memnum = cell(1,2);
for i=1:2
    memnum{i} = mean(squeeze(result.output.memnumbytime(1:2:end-1,n*i,3)));
end
save(outputFile, 'memnum');
end
