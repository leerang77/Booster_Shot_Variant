fileList = dir(fullfile('..','parameters'));
a = arrayfun(@(x) x.name, fileList, 'UniformOutput', false);
fid = fopen('jobs_to_run.txt', 'w');
for i=1:length(a)
    if a{i}(1)=='V'
        fprintf(fid, ['./RunBatchJob < ../parameters/',a{i},'\n']);
    end
end
fclose(fid);