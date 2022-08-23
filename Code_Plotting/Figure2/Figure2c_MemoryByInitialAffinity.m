fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%
[result, param] = loadFigure2Data();
result{2} = result{3}; param{2} = param{3};
num_by_E0_sum = cell(1,2);
for i=1:2
    num_by_E0_sum{i} = numByInitialAff(result{i});
end
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'num_by_E0_sum')

%% Number by initial affinity

function num_by_E0_sum = numByInitialAff(result)
f = 7;
output = result.output.finalmem;
naiveaff = squeeze(result.naive(:,:,3));
n_GC = size(output,1);
aff = 6:0.2:8.2;
num_by_E0 = zeros(2, n_GC, length(aff)-1);
for k=1:n_GC
   affs = squeeze(output(k,:,3));
   lineages = squeeze(output(k,affs>f,1));
   targets = squeeze(output(k,affs>f,2));
   E0 = naiveaff(k,lineages);
   for target=1:2
       num_by_E0(target,k,:) = histcounts(E0(targets==target), aff);
   end
end
num_by_E0_sum = squeeze(sum(num_by_E0, 2));
sum(sum(num_by_E0_sum(:,6:end)))/sum(sum(num_by_E0_sum))
num_by_E0_sum(num_by_E0_sum==0)=nan;
% xticklabels(cellfun(@(x) num2str(x, '%.1f'), num2cell(6:0.2:8), 'UniformOutput', false))
end