fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
p = base_case_parameters();
result = cell(1,2); param = cell(1,2);
[result{1},dirname] = load_result(p);
param{1} = result{1}.param;

p.vaxnum = 2; p.tmax = 180;
[result{2},dirname] = load_result(p);
param{2} = result{2}.param;

%%
naive_entry = repmat({zeros(2,22)},1,2);
for ep=1:2
for t=0:21
    for vaxcase=1:2
        naive_entry{ep}(vaxcase,1+t) = mean(sum(squeeze(result{vaxcase}.naive(:,:,5)>0 &...
         result{vaxcase}.naive(:,:,5)<t & result{vaxcase}.naive(:,:,2)==ep),2) );
    end
end
end

t = 8;
affs = 6:0.2:8;
num_by_aff = zeros(2,2,length(affs));
for ep=1:2
for vaxcase = 1:2
    idx = result{vaxcase}.naive(:,:,5)>0 & result{vaxcase}.naive(:,:,5)<t...
        & result{vaxcase}.naive(:,:,2)==ep;
    naive_affs = result{vaxcase}.naive(:,:,3);
    naive_entry_affs = squeeze(naive_affs(idx));
    num_by_aff(ep,vaxcase,:) = histcounts(naive_entry_affs, [affs,inf]);
end
end

num_by_aff = num_by_aff/2000;
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile, 'affs', 'num_by_aff', 'naive_entry', 'param')


