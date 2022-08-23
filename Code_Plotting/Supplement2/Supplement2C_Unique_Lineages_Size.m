fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
% GC
p = base_case_parameters();
p.last = 200;
result = cell(1,2);
[result{1},~] = load_result(p);
p.vaxnum = 2; p.tmax = 180;
[result{2},~] = load_result(p);
GC = cell(1,2); EGC = cell(1,2);
f = 7;
for i=1:2
    output = result{i}.output.finalmem;
    n_GC = size(output,1);
    lin = 1:2011;
    num_by_lineage = zeros(n_GC, 2010);
    for k=1:n_GC
       affs = squeeze(output(k,:,3));
       lineages = squeeze(output(k,affs>f,1));
       num_by_lineage(k,:) = histcounts(lineages, lin);
    end
    lineages_size = reshape(num_by_lineage,1,[]);
    GC{i} = sort(lineages_size(lineages_size>0),'descend');
end

% EGC
uniqueidx = result{2}.finalmemEGC(9,:);
uniqueidx = uniqueidx(result{2}.finalmemEGC(3,:)>f);
cnts = arrayfun(@(x) sum(x==uniqueidx), unique(uniqueidx));
EGC{2} = sort(cnts, 'descend');


outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'GC', 'EGC')


